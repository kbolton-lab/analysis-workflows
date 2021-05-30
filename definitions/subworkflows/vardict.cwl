#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "vardict parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
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
        type: string?
        default: 'TUMOR'
    normal_sample_name:
        type: string
        default: 'NORMAL'
    af_threshold:
        type: float?
        default: 0.005
outputs:
    unfiltered_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    # unfiltered_vcf:
    #     type: File
    #     outputSource: filter/unfiltered_vcf
    #     secondaryFiles: [.tbi]
    # filtered_vcf:
    #     type: File
    #     outputSource: filter/filtered_vcf
    #     secondaryFiles: [.tbi]
steps:
    split_interval_list:
        run: ../tools/split_interval_list_to_bed.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_beds]
    vardict:
        scatter: interval_list
        run: ../tools/vardict.cwl
        in:
            reference: reference
            af_threshold: af_threshold
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list/split_beds
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
            output_name: 
                valueFrom: $(inputs.tumor_sample_name).vardict.vcf
        out:
            [vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: vardict/vcf
        out:
            [merged_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
    # filter:
    #     run: fp_filter.cwl
    #     in:
    #         reference: reference
    #         bam: tumor_bam
    #         vcf: index/indexed_vcf
    #         variant_caller: 
    #             valueFrom: "vardict"
    #         sample_name: tumor_sample_name
    #     out:
    #         [unfiltered_vcf, filtered_vcf]
    # bcbio_filter:
    #     run: ../tools/bcftools_filter.cwl



