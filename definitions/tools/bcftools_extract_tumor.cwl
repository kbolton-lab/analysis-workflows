#!/usr/bin/env cwl-runner

### Need to do general isec ###
cwlVersion: v1.0
class: CommandLineTool
label: "extract tumor sample from tumor/normal"

baseCommand: ["/opt/bcftools/bin/bcftools", "view"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: InlineJavascriptRequirement

inputs:
    vcf:
        type: File
        inputBinding:
            position: 2
        doc: "input bgzipped tabix indexed vcf to obtain complement"
        secondaryFiles: [.tbi]
    tumor_sample_name:
        type: string
        inputBinding:
            position: 1
            prefix: "-s"
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 3
            prefix: "--output-type"
        doc: "output file format"
    output_vcf_name:
        type: string?
        default: "tumor.vcf.gz"
        inputBinding:
            position: 4
            prefix: "--output-file"
        doc: "output vcf file name"
    
outputs:
    tumor_only_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
