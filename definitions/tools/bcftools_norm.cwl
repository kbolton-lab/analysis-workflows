#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/opt/bcftools/bin/bcftools", "norm"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"


inputs:
    multiallelic:
        type: string?
        default: "-any"
        inputBinding:
            position: 1
            prefix: "--multiallelics"
        doc: "rsplit multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+) [snps|indels|both|any]"
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
        default: "bcftools_norm.vcf.gz"
        inputBinding:
            position: 4
            prefix: "--output"
        doc: "output vcf file name"
    vcf:
        type: File
        inputBinding:
            position: 2
        doc: "input bgzipped tabix indexed vcf to normalize"

outputs:
    normalized_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

