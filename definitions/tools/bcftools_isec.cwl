#!/usr/bin/env cwl-runner

### Need to do ###
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/opt/bcftools/bin/bcftools", "norm"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"

inputs:
    force_merge:
        type: boolean?
        default: true
        inputBinding:
            position: 1
            prefix: "--force-samples"
        doc: "resolve duplicate sample names"
    merge_method:
        type:
            type: string
        default: "-any"
        inputBinding:
            position: 2
            prefix: "--multiallelics"
        doc: "split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). [snps|indels|both|any]"
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
            position: 5
        doc: "input bgzipped tabix indexed vcf to normalize"

outputs:
    norm_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

