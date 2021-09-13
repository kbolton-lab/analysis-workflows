#!/usr/bin/env cwl-runner

### Need to do general isec ###
cwlVersion: v1.0
class: CommandLineTool
label: "Used for obtaining variants by intersecting vcfs, returns first vcf perspective"

baseCommand: ["/opt/bcftools/bin/bcftools", "isec"]
arguments: [
    { position: 1, valueFrom: "-n+2" },
    { position: 2, valueFrom: "-w1" },
]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"

inputs:
    vcf1:
        type: File
        inputBinding:
            position: 3
        doc: "input bgzipped tabix indexed vcf to obtain complement"
        secondaryFiles: [.tbi]
    vcf2:
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
    overlap_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
