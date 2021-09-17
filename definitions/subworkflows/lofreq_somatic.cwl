#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "lofreq somatic workflow"
inputs:
    reference:
        type: File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai]
    roi_bed:
        type: File?
    interval_list:
        type: File
    scatter_count:
        type: int
    tumor_sample_name:
        type: string
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
    split_interval_list_to_bed:
        run: ../tools/split_interval_list_to_bed.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_beds]
    lofreq_somatic:
        scatter: interval_list
        run: ../tools/lofreq_somatic.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list_to_bed/split_beds
        out:
            [vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: lofreq_somatic/vcf
        out:
            [merged_vcf]
    reformat_vcf:
        run: ../tools/lofreq_reformat.cwl
        in:
            vcf: merge/merged_vcf
            tumor_sample_name: tumor_sample_name
        out:
            [reformat_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: reformat_vcf/reformat_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: index/indexed_vcf
            variant_caller:
                valueFrom: "lofreq"
            sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]