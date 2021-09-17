#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "ArcherDX FastQ Umi Extraction"
requirements:
    - class: SchemaDefRequirement
      types:
        - $import: ../types/sequence_data.yml
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: InlineJavascriptRequirement
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data
    umi_length:
        type: int
        default: 8
outputs:
    fastqs:
        type: File[]
        outputSource: repair/fastqs
    fastq1:
        type: File
        outputSource: repair/fastq1
    fastq2:
        type: File
        outputSource: repair/fastq2
steps:
    filter_umi_length:
        run: ../tools/filter_umi_length.cwl
        in:
            fastq1:
                source: sequence
                valueFrom: "$(self.sequence.hasOwnProperty('fastq1') self.sequence.fastq1 : null)"
            fastq2:
                source: sequence
                valueFrom: "$(self.sequence.hasOwnProperty('fastq2') self.sequence.fastq2 : null)"
            umi_length: umi_length
        out:
            [fastq1, fastq2]
    repair:
        run: ../tools/bbmap_repair.cwl
        in:
            fastq1:
                source: filter_umi_length/fastq1
            fastq2:
                source: filter_umi_length/fastq2
        out:
            [fastqs, fastq1, fastq2]
