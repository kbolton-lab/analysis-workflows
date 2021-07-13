#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Used for complementing first vcf from second vcf"

baseCommand: ["/usr/local/bin/bcftools", "isec"]
arguments: [
    { position: 1, valueFrom: "-C" },
    { position: 2, valueFrom: "-w1" },
]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"

inputs:
    vcf:
        type: File
        inputBinding:
            position: 3
        doc: "input bgzipped tabix indexed vcf to obtain complement"
        secondaryFiles: [.tbi]
    exclude_vcf:
        type: File
        inputBinding:
            position: 4
        doc: "input bgzipped tabix indexed vcf to exclude"
        secondaryFiles: [.tbi]
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 5
            prefix: "--output-type"
        doc: "output file format"
    output_vcf_name:
        type: string?
        default: "bcftools_isec.vcf.gz"
        inputBinding:
            position: 6
            prefix: "--output"
        doc: "output vcf file name"
    

outputs:
    complement_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

