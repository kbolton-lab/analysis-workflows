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
    gnomad_AF_only:
        type: File
        secondaryFiles: [.tbi]
    gnomad_AF:
        type: float
        default: .005
    filter_flag:
        type:
          - "null" 
          - type: enum
            symbols: ["include", "exclude"]
        default: "include"
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    normal_bams:
        type: 
            type: array
            items: File
        secondaryFiles: [.bai]
    mapq:
        type: int?
        default: 0
    baseq:
        type: int?
        default: 0
    pon_final_name:
        type: string?
        default: "pon.total.counts"
        doc: "prefix final name for total pileup, usually caller name + 'pon.total.counts'"
    pon_pvalue:
        type: string?
        default: "1.293211e-09"
outputs:
    processed_gnomAD_vcf:
        type: File
        outputSource: index_pon_vcf/indexed_vcf
        doc: "only gnomad filter"
        secondaryFiles: [.tbi]
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
    gnomAD_filter:
        run: ../subworkflows/gnomAD_filter.cwl
        in:
            vcf: caller_vcf
            gnomad_AF_only: gnomad_AF_only
            gnomad_AF: gnomad_AF
            filter_flag: filter_flag
            output_prefix: caller_prefix
        out:
            [filtered_vcf]
        doc: "hard filter, output vcf format not .gz for pon"     
    pon_pileup:
        run: get_base_counts_multi_sample.cwl
        in:
            reference: reference
            normal_bams: normal_bams
            vcf: gnomAD_filter/filtered_vcf
            mapq: mapq
            baseq: baseq
            pon_final_name: pon_final_name
        out:
            [total_counts]
    call_R_fisher: 
        run: ../tools/normal_fisher.cwl
        in:
            vcf: gnomAD_filter/filtered_vcf
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
    