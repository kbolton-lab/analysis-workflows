#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "mutect parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: [^.bai, .bai]
    normal_bam:
        type: File?
        secondaryFiles: [^.bai]
    interval_list:
        type: File
    scatter_count:
        type: int
    tumor_sample_name:
        type: string
    min_var_freq:
        type: float?
        default: 0.05
outputs:
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
steps:
    mutect:
        run: ../tools/mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
        out:
            [vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: mutect/vcf
            variant_caller:
                valueFrom: "mutect"
            sample_name: tumor_sample_name
            min_var_freq: min_var_freq
        out:
            [unfiltered_vcf, filtered_vcf]
