#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Using pyvcf count number of samples with variant in a multivcf. https://pyvcf.readthedocs.io/en/latest/FILTERS.html"
doc: "If the number of samples with variant is greater than (>) these are flagged and then removed"

baseCommand: ["/bin/bash", "PON.sh"]
arguments: [
    { position: 5, valueFrom: $(inputs.vcf.nameroot) }
]

requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'PON.sh'
        entry: |
            
            set -eou pipefail
            
            export caller="$3"
            export vcf_in="$1"
            export pon_total="$2"
            export name="$5"
            export p_value="$4"

            printf "##INFO=<ID=PON_RefDepth,Number=1,Type=Integer,Description=\"Total Ref_Depth for Normals\">\n##INFO=<ID=PON_AltCounts,Number=1,Type=Integer,Description=\"Total Alt_Counts for Normals\">\n" > pileup.header;
            printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
            printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name \(with whitespace translated to underscores)\">" > sample.header;

            sample=`bcftools query -l $vcf_in`
            bcftools view -H $vcf_in | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
            bgzip $sample.name;
            tabix $sample.name.gz -s1 -b2 -e2;
            bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf_in -Ov -o $name.sample.vcf;
            
            ## the filter_nsamples.py python script has to be in the current directory
            if [[ $caller =~ "[Vv]arscan"  ]]
            then
                echo '
                #!/usr/bin/env Rscript

                args = commandArgs(trailingOnly=TRUE)

                if (length(args)==0) {
                stop("At least one argument must be supplied (input file).n", call.=FALSE)
                } else if (length(args)==1) {
                # default output file
                stop("Must supply (output file).n", call.=FALSE)
                }

                df = read.table(args[1], header=F)
                
                if (length(colnames(df)) != 8) {
                stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_AltCounts\t%INFO/PON_RefDepth\t[%AD]\t[%RD]", call.=FALSE)
                }

                df$fisher.exact.pval <- apply(df, 1, function(x) {
                x <- as.numeric(x[-c(1,2,3,4)])
                if ((x[1]+x[2])==0) {
                    return(0)
                } else {
                    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
                }
                })
                write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
                ' > fisherTestInput.R
                bcftools annotate -a $pon_total -h pileup.header -c CHROM,POS,REF,ALT,PON_RefDepth,PON_AltCounts $name.sample.vcf -Ov -o $name.sample.pileup.vcf;
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_AltCounts\t%INFO/PON_RefDepth\t[%AD]\t[%RD]\n' $name.sample.pileup.vcf > $name.fisher.input;
            else
                echo '
                #!/usr/bin/env Rscript 

                args = commandArgs(trailingOnly=TRUE)

                if (length(args)==0) {
                stop("At least one argument must be supplied (input file).n", call.=FALSE)
                } else if (length(args)==1) {
                # default output file
                stop("Must supply (output file).n", call.=FALSE)
                }

                df = read.table(args[1], header=F)
                
                #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
                df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
                if (length(colnames(df)) != 8) {
                stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_AltCounts\t%INFO/PON_RefDepth\t[%AD]", call.=FALSE)
                }

                df$fisher.exact.pval <- apply(df, 1, function(x) {
                x <- as.numeric(x[-c(1,2,3,4)])
                if ((x[1]+x[2])==0) {
                    return(0)
                } else {
                    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
                }
                })
                write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
                ' > fisherTestInput.R
                 bcftools annotate -a $pon_total -h pileup.header -c CHROM,POS,REF,ALT,PON_RefDepth,PON_AltCounts $name.sample.vcf -Ov -o $name.sample.pileup.vcf;
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltCounts\t[%AD]\n' $name.sample.pileup.vcf > $name.fisher.input;
                
            fi
            chmod u+x fisherTestInput.R

            
            LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
            bgzip -f $name.fisher.output
            tabix -f -s1 -b2 -e2 $name.fisher.output.gz
            bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $name.sample.pileup.vcf -Ov -o $name.pileup.fisherPON.vcf
            bcftools filter -i "INFO/PON_FISHER<$p_value" $name.pileup.fisherPON.vcf -Ov -o $name.filtered.pileup.fisherPON.vcf

inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    pon_total:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
    caller:
        type: string
        default: "mutect2"
        inputBinding:
            position: 3
    p_value:
        type: string?
        default: ".05"
        inputBinding:
            position: 4
    
outputs:
    pon_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.vcf.nameroot).pileup.fisherPON.vcf"
    pon_filtered_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.vcf.nameroot).filtered.pileup.fisherPON.vcf"


