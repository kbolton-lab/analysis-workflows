#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "mutect parallel workflow"
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
        secondaryFiles: [^.bai, .bai]
    interval_list:
        type: File
    scatter_count:
        type: int
    tumor_sample_name:
        type: string
outputs:
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
    unfiltered_vcf:
        type: File
        outputSource: norm_index/indexed_vcf
        secondaryFiles: [.tbi]
        doc: "This is the unfiltered from fp_filter.cwl that is normalized with bcftools for use in nsamples.cwl"
steps:
    split_interval_list:
        run: ../tools/split_interval_list.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_interval_lists]
    mutect:
        scatter: interval_list
        run: ../tools/mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list/split_interval_lists
        out:
            [vcf, tumor_sample_name]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: mutect/vcf
        out:
            [merged_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: index/indexed_vcf
            variant_caller: 
                valueFrom: "mutect"
            sample_name: tumor_sample_name
            # sample_name: 
            #     source: mutect/tumor_sample_name
            #     valueFrom: |
            #         ${
            #             return self[0]
            #         }
        out:
            [unfiltered_vcf, filtered_vcf]
    bcftools_norm:
        run: ../tools/bcftools_norm.cwl
        in: 
            vcf: filter/unfiltered_vcf
        out:
            [normalized_vcf]
        doc: "required so that we can split multi-allelic in a better way than GATK does for counting Nsamples by pyvcf"
    norm_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bcftools_norm/normalized_vcf
        out:
            [indexed_vcf]

    