#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Run gnomad filter and PoN pileup filter"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    caller_vcf:
        type: File
        secondaryFiles: [.tbi]
    caller_prefix: 
        type: string?
        default: "caller"
    gnomAD_exclude_vcf:
        type: File
        secondaryFiles: [.tbi]
    normal_bams:
        type: 
            type: array
            items: string
    mapq:
        type: int?
        default: 5
    baseq:
        type: int?
        default: 5
    pon_final_name:
        type: string?
        default: "pon.total.counts"
        doc: "prefix final name for total pileup, usually caller name + 'pon.total.counts'"
    pon_pvalue:
        type: string?
        default: "1.282123e-09"
outputs:
    processed_filtered_vcf:
        type: File
        outputSource: index_pon_filtered_vcf/indexed_vcf
        doc: "gnomad and pon filter"
        secondaryFiles: [.tbi]
    pon_total_counts:
        type: File
        outputSource: pon_pileup/total_counts
        secondaryFiles: [.tbi]

steps: 
    pon_pileup:
        run: get_base_counts_multi_sample.cwl
        in:
            reference: reference
            normal_bams: normal_bams
            vcf: isec_complement_gnomAD/complement_vcf
            mapq: mapq
            baseq: baseq
            pon_final_name: pon_final_name
        out:
            [total_counts]
    call_R_fisher: 
        run: ../tools/normal_fisher.cwl
        in:
            vcf: isec_complement_gnomAD/complement_vcf
            pon_total: pon_pileup/total_counts
            p_value: pon_pvalue
            caller: caller_prefix
        out:
            [pon_vcf, pon_filtered_vcf]
        doc: "hard and soft filter but only annotating hard filter with VEP"
    ##add sample name as INFO tag
    index_pon_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: call_R_fisher/pon_vcf
        out:
            [indexed_vcf]
    index_pon_filtered_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: call_R_fisher/pon_filtered_vcf
        out:
            [indexed_vcf]   
    