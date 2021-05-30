#!/usr/bin/env cwl-runner

### Need to do general isec ###
cwlVersion: v1.0
class: CommandLineTool
label: "Used for complementing first vcf from second vcf"

baseCommand: ["/opt/bcftools/bin/bcftools", "view"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"

inputs:
    vcf:
        type: File
        inputBinding:
            position: 3
        doc: "input bgzipped tabix indexed vcf to obtain complement"
        secondaryFiles: [.tbi]
    tumor_sample_name:
        type: string
        inputBinding:
            position: 2
            prefix: "-s"
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 4
            prefix: "--output-type"
        doc: "output file format"
    
outputs:
    tumor_only_vcf:
        type: File
        outputBinding:
            glob: $(inputs.tumor_sample_name).vcf.gz
