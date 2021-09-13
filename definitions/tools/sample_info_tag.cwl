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
      ramMin: 2000
      tmpdirMin: 100000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sample.sh'
        entry: |
            
            set -eou pipefail
            
        
            export vcf_in="$1"
            export output_name="$2"
            export output_type="$3"
           
            printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

            sample=`bcftools query -l $vcf_in`
            bcftools view -H $vcf_in | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
            bgzip $sample.name;
            tabix $sample.name.gz -s1 -b2 -e2;
            bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf_in $output_type $output_name;
            tabix $output
            
            

inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    output_vcf_name:
        type: string
        default: "sample_info.vcf.gz"
        inputBinding:
            position: 2
        doc: "output vcf file name"
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 3
            prefix: "--output-type"
        doc: "output file format"
    
outputs:
    sample_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)


