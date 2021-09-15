#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "lofreq v2.1.3.1 call"
baseCommand: ["/bin/bash", "lofreq_reformat_helper.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 2000
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'lofreq_reformat_helper.sh'
        entry: |
            #!/bin/bash

            set -o errexit
            set -o nounset

            if [ $# -lt 2 ]
            then
                echo "Usage: $0 [VCF] [SAMPLE_NAME]"
                exit 1
            fi

            zcat $1 | grep "##" > lofreq.reformat.vcf;
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$2" >> lofreq.reformat.vcf;
            zcat $1 | grep -v '#' | awk '{ n=split($8, semi, /;/); sample=""; format=""; for(i in semi){ split(semi[i], equ, /=/); if(i<=3){ if(i+1==4) sample=sample equ[2]; else sample=sample equ[2] ":"; if(i+1==4) format=format equ[1]; else format=format equ[1] ":";}}{print $0, format, sample}}' OFS='\t' >> lofreq.reformat.vcf;
            bgzip lofreq.reformat.vcf;
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    tumor_sample_name:
        type: string
        inputBinding:
            position: 2
outputs:
    reformat_vcf:
        type: File
        outputBinding:
            glob: "lofreq.reformat.vcf.gz"
