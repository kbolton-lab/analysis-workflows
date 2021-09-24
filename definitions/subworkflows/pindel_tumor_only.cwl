#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "pindel parallel workflow"
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
        secondaryFiles: [.bai, ^.bai]
    interval_list:
        type: File
    insert_size:
        type: int
        default: 400
    scatter_count:
        type: int
        default: 50
    tumor_sample_name:
        type: string
    ref_name:
        type: string?
        default: "GRCh38DH"
    ref_date:
        type: string?
        default: "20161216"
    pindel_min_supporting_reads:
        type: int?
        default: 3
outputs:
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [.tbi]
        doc: "This is the unfiltered from fp_filter.cwl that is normalized with bcftools for use in nsamples.cwl"
steps:
    split_interval_list_to_bed:
        run: ../tools/split_interval_list_to_bed.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out:
            [split_beds]
    pindel_cat:
        scatter: region_file
        run: pindel_cat.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            region_file: split_interval_list_to_bed/split_beds
            insert_size: insert_size
            tumor_sample_name: tumor_sample_name
        out:
            [per_region_pindel_out]
    cat_all:
        run: ../tools/cat_all.cwl
        in:
            region_pindel_outs: pindel_cat/per_region_pindel_out
        out:
            [all_region_pindel_head]
    pindel2vcf:
        run: ../tools/pindel2vcf.cwl
        in:
            reference: reference
            pindel_out: cat_all/all_region_pindel_head
            ref_name: ref_name
            ref_date: ref_date
            min_supporting_reads: pindel_min_supporting_reads
        out:
            [pindel_vcf]
    bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: pindel2vcf/pindel_vcf
        out:
            [bgzipped_file]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
    remove_end_tags:
        run: ../tools/remove_end_tags.cwl
        in:
            vcf: index/indexed_vcf
        out:
            [processed_vcf]
    reindex:
        run: ../tools/index_vcf.cwl
        in:
            vcf: remove_end_tags/processed_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: reindex/indexed_vcf
            variant_caller:
                valueFrom: "pindel"
            sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
