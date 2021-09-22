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
    interval_list:
        type: File
    scatter_count:
        type: int
    tumor_sample_name:
        type: string?
        default: 'TUMOR'
    af_threshold:
        type: float?
        default: 0.005
    bcbio_filter_string:
        type: string
        default: "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 45)))"
        doc: "http://bcb.io/2016/04/04/vardict-filtering/"
outputs:
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [.tbi]
    bcbio_filtered_vcf:
        type: File
        #outputSource: filter/filtered_vcf
        outputSource: index_bcbio/indexed_vcf
        secondaryFiles: [.tbi]
        doc: "this is the bcbio filtered from the fp_filter filter/unfiltered_vcf"
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
            interval_list: split_interval_list/split_beds
            tumor_sample_name: tumor_sample_name
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
    filter:
        run: fp_filter_vardict.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: index/indexed_vcf
            variant_caller:
                valueFrom: "vardict"
            sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    bcftools_norm:
        run: ../tools/bcftools_norm.cwl
        in:
            reference: reference
            vcf: filter/unfiltered_vcf
            output_vcf_name:
                valueFrom: "vardict_full.vcf.gz"
        out:
            [normalized_vcf]
        doc: "required so that we can split multi-allelic in a better way than GATK does for counting Nsamples by pyvcf"
    norm_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bcftools_norm/normalized_vcf
        out:
            [indexed_vcf]
    bcbio_filter:
        run: ../tools/bcftools_filter_bcbio.cwl
        in:
            vcf: norm_index/indexed_vcf
            filter_string: bcbio_filter_string
            filter_flag:
                valueFrom: "exclude"
            output_type:
                default: "z"
            output_vcf_name:
                valueFrom: $(inputs.vcf.nameroot.replace(".vcf","")).bcbiofilter.vcf.gz
        out: [filtered_vcf]
        doc: "this is bcbio's filter for vardict"
    index_bcbio:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bcbio_filter/filtered_vcf
        out:
            [indexed_vcf]
