#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "lofreq somatic workflow"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
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
    lofreq_vcf:
        type: File
        outputSource: lofreq/vcf
        secondaryFiles: [.tbi]
steps:
    split_interval_list_to_bed:
        run: ../tools/split_interval_list_to_bed.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_beds]
    lofreq:
        scatter: interval_list
        run: ../tools/lofreq_call.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            roi_bed: split_interval_list_to_bed/split_beds
        out:
            [vcf]
