#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "add_sample_name_tag.sh"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "sam_bcftools_tabix_bgzip:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "add_sample_name_tag.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          basen=`basename "$1"`
          basen="$2.sample.$basen"

          printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header
          sample=$(bcftools query -l "$1")

          /usr/local/bin/bcftools view -H "$1" | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
          /usr/local/bin/bgzip $sample.name;
          /usr/local/bin/tabix $sample.name.gz -s1 -b2 -e2;
          /usr/local/bin/bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLENAME_TAG "$1" -o "$basen";

inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 1
        doc: "vcf file to filter"
    caller:
        type: string
        inputBinding:
            position: 2
        doc: "Name of caller to prepend vcf file"

outputs:
    renamed_vcf:
        type: File