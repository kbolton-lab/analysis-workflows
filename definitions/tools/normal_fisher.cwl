#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run a fisher test on the PoN bam totals for each variant"


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


            printf "##INFO=<ID=PON_RefDepth,Number=1,Type=Integer,Description=\"Total Ref_Depth for Normals\">\n##INFO=<ID=PON_AltDepth,Number=1,Type=Integer,Description=\"Total Alt_Counts for Normals\">\n##INFO=<ID=PON_VAF,Number=1,Type=Float,Description=\"VAF for Normals\">\n" > pileup.header;
            printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
            printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

            sample=`bcftools query -l $vcf_in`
            bcftools view -H $vcf_in | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
            bgzip $sample.name;
            tabix $sample.name.gz -s1 -b2 -e2;
            bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf_in -Ov -o $name.sample.vcf;

            ## Varscan has AD and RD instead of comma sep AD field
            ## you can't double quote for string match in bash like you can in zsh so need to make it a variable
            pat="[Vv]arscan"

            ## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
            patt="[Ll]ofreq"
            if [[ $caller =~ $pat ]]
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
                stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]", call.=FALSE)
                }

                df$fisher.exact.pval <- apply(df, 1, function(x) {
                x <- as.numeric(x[-c(1,2,3,4)])
                if (x[1]+x[2]==0){
                    return(0)
                } else if (x[2]==0 & x[1]!=0) {
                    return(0)
                } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                    return(1)
                } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+ x[4])) {
                    return(1)
                } else {
                    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
                }
                })
                write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
                ' > fisherTestInput.R
                bcftools annotate -a $pon_total -h pileup.header -c CHROM,POS,REF,ALT,PON_RefDepth,PON_AltDepth,PON_VAF $name.sample.vcf -Ov -o $name.sample.pileup.vcf;
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]\n' $name.sample.pileup.vcf > $name.fisher.input;
            elif [[ $caller =~ $patt ]]
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

                #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
                df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]

                if (length(colnames(df)) < 8) {
                stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4", call.=FALSE)
                }

                df$fisher.exact.pval <- apply(df, 1, function(x) {
                x <- as.numeric(x[-c(1,2,3,4)])
                # Remember, Lofreq splits DP4 into RefFwd, RefRev and AltFwd, AltRev so technically ref = x[3] + x[4] and alt = x[5] + x[6]
                ref = x[3] + x[4]
                alt = x[5] + x[6]
                if (x[1]+x[2]==0){
                    return(0)
                } else if (x[2]==0 & x[1]!=0) {
                    return(0)
                } else if ((x[1]==0 & x[2]!=0) | (ref==0 & alt!=0)) {
                    return(1)
                } else if (x[2]/(x[1]+x[2]) >= alt/(ref+alt)) {
                    return(1)
                } else {
                    return(fisher.test(matrix(c(x[1], x[2], ref, alt), ncol=2))$p.value)
                }
                })
                write.table(df[, -c(9:10)], file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
                ' > fisherTestInput.R
                bcftools annotate -a $pon_total -h pileup.header -c CHROM,POS,REF,ALT,PON_RefDepth,PON_AltDepth,PON_VAF $name.sample.vcf -Ov -o $name.sample.pileup.vcf;
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf > $name.fisher.input;
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
                stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]", call.=FALSE)
                }

                df$fisher.exact.pval <- apply(df, 1, function(x) {
                x <- as.numeric(x[-c(1,2,3,4)])
                if (x[1]+x[2]==0){
                    return(0)
                } else if (x[2]==0 & x[1]!=0) {
                    return(0)
                } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                    return(1)
                } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
                    return(1)
                } else {
                    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
                }
                })
                write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
                ' > fisherTestInput.R
                bcftools annotate -a $pon_total -h pileup.header -c CHROM,POS,REF,ALT,PON_RefDepth,PON_AltDepth,PON_VAF $name.sample.vcf -Ov -o $name.sample.pileup.vcf;
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $name.sample.pileup.vcf > $name.fisher.input;

            fi
            chmod u+x fisherTestInput.R


            LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
            bgzip -f $name.fisher.output
            tabix -f -s1 -b2 -e2 $name.fisher.output.gz
            bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $name.sample.pileup.vcf -Oz -o $name.pileup.fisherPON.vcf.gz
            bcftools filter -i "INFO/PON_FISHER<$p_value" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz

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
            glob: "$(inputs.vcf.nameroot).pileup.fisherPON.vcf.gz"
    pon_filtered_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.vcf.nameroot).filtered.pileup.fisherPON.vcf.gz"
