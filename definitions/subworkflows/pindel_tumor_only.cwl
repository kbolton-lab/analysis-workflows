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
    somaticfilter:
        run: ../tools/pindel_somatic_filter.cwl
        in:
            reference: reference
            pindel_output_summary: cat_all/all_region_pindel_head
        out:
            [vcf]
    bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: somaticfilter/vcf
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
    bcftools_norm:
        run: ../tools/bcftools_norm.cwl
        in:
            reference: reference
            vcf: filter/unfiltered_vcf
            output_vcf_name:
                valueFrom: "pindel_full.vcf.gz"
        out:
            [normalized_vcf]
        doc: "required so that we can split multi-allelic in a better way than GATK does for counting Nsamples by pyvcf"
    norm_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bcftools_norm/normalized_vcf
        out:
            [indexed_vcf]
