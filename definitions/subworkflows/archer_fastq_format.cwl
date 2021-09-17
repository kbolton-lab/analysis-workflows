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
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
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
        scatter: sequence
        run: ../tools/filter_umi_length.cwl
        in:
            sequence: sequence
            umi_length: umi_length
        out:
            [fastq1, fastq2]
    repair:
        run: ../tools/bbmap_repair.cwl
        in:
            fastq1:
                source: filter_umi_length/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: filter_umi_length/fastq2
                linkMerge: merge_flattened
        out:
            [fastqs, fastq1, fastq2]
