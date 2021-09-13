#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/opt/bcftools/bin/bcftools", "filter"]

arguments: [
    { position: 6, valueFrom: "BCBIO", prefix: "-s" },
    "-m+"
]

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
            position: 3
        doc: "input bgzipped tabix indexed vcf to filter"
    filter_flag:
        type:
            type: enum
            symbols: ["include", "exclude"]
        default: "include"
        inputBinding:
            position: 1
            valueFrom: |
                ${
                    if(inputs.filter_flag === "include") {
                        return "-i";
                    } else {
                        return "-e";
                    }
                }    
    filter_string:
        type: string
        inputBinding:
            position: 2         
    output_vcf_name:
        type: string?
        default: "bcftools_filter.vcf.gz"
        inputBinding:
            position: 5
            prefix: "--output"
        doc: "output vcf file name"
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
    filtered_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

