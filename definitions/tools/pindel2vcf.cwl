#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
baseCommand: ["/usr/bin/pindel2vcf"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.3.1"
arguments:
    ["-G"]
inputs:
    pindel_out:
        type: File
        inputBinding:
            prefix: "-p"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-r"
    ref_name:
        type: string?
        inputBinding:
            prefix: "-R"
        default: "GRCh38DH"
    min_supporting_reads:
        type: int?
        inputBinding:
            prefix: "-e"
        default: 3
    ref_date:
        type: string?
        inputBinding:
            prefix: "-d"
        default: "20161216"
    output_name:
        type: string?
        inputBinding:
            prefix: "-v"
        default: "pindel.vcf"
outputs:
    pindel_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.output_name)"
