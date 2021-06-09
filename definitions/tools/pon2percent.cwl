#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "PON2.sh"]
arguments: [
    { position: 4, valueFrom: $(inputs.caller).final.annotated.vcf.gz }
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
      - entryname: 'PON2.sh'
        entry: |
            
            set -eou pipefail
            
            export caller="$3"
            export vcf_in="$1"
            export normal2="$2"
            export name="$4"

            printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
            bcftools query -f "%CHROM\t%POS\t%REF\t%REF\t1\n" $normal2 > normal2.txt
            bgzip -f normal2.txt
            tabix -f -s1 -b2 -e2 normal2.txt.gz
            bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent $vcf_in -Oz -o $name
            tabix $name


inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    vcf2PON:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
        doc: "this is variant calls from PoN with same variant caller as tumor"
    caller:
        type: string
        default: "mutect2"
        inputBinding:
            position: 3
    
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.caller).final.annotated.vcf.gz"
        secondaryFiles: [.tbi]


