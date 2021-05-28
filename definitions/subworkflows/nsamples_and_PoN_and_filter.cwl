#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Run NSamples and Pileup creation"
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
    caller_vcfs:
        type: 
            type: array
            items: File
        secondaryFiles: [.tbi]
    caller_prefix: 
        type: string?
        default: "caller"
    merge_method:
        type: 
          - "null"
          - type: enum
            symbols: ["none", "snps", "indels", "both", "all", "id"]
        default: "none"
    gnomad_AF_only:
        type: File
        secondaryFiles: [.tbi]
    gnomad_AF:
        type: float
        default: .005
    filter_flag:
        type:
            type:
            - "null" 
            - type: enum
              symbols: ["include", "exclude"]
            default: "include"
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    nsamples:
        type: int
        default: 5
        doc: |
            The threshold for number of samples a mutation can have to filter out artefacts and sequencing errors
            Can be found with qbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE) in R. If number of samples with variant > nsamples these are removed
    normal_bams:
        type: 
            type: array
            items: File
        secondaryFiles: [.bai]
        type: int
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
outputs:
    final_concatenated_vcf:
        type: File
        outputSource: concat_vcfs/pon_concatenated_vcf
    final_concatenated_filtered_vcf:
        type: File
        outputSource: concat_vcfs/pon_filtered_concatenated_vcf
    merged_nsamples:
        type: File
        outputSource: nsamples/merged_nsamples
    merged_nsamples_underN:
        type: File
        outputSource: nsamples/merged_nsamples_underN
    merged_nsamples_overN:
        type: File
        outputSource: nsamples/merged_nsamples_overN
        secondaryFiles: [.tbi]
    pon_total_counts:
        type: File
        outputSource: pon_pileup/total_counts
        secondaryFiles: [.tbi]

steps:
    merge:
        run: ../tools/bcftools_merge.cwl
        in:
            vcfs: caller_vcfs
            merge_method: merge_method
            output_vcf_name: 
                source: "#caller_prefix"
                valueFrom: "$(self).merged.vcf.gz"
            merge_method: merge_method
            output_type: "z"
        out:
            [merged_sv_vcf]
    merge_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge/merged_sv_vcf
        out:
            [indexed_vcf]
    gnomAD_filter:
        run: ../subworkflows/gnomAD_filter.cwl
        in:
            vcf: merge_index/indexed_vcf
            gnomad_AF_only: gnomad_AF_only
            gnomad_AF: gnomad_AF
            flag: filter_flag
            output_prefix: caller_prefix
        out:
            [vcf_below_gnomad_AF]
        doc: "hard filter"
    gnomAD_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: gnomAD_filter/vcf_below_gnomad_AF
        out:
            [indexed_vcf]       
    nsamples: # done
        run: ../tools/nsamples.cwl
        in:
            output_name: 
                source: "#caller_prefix"
                valueFrom: "$(self).merged.nsamples"
            merged_vcf: gnomAD_index/indexed_vcf
            nsamples: nsamples
        out:
            [merged_nsamples, merged_nsamples_underN, merged_nsamples_overN]
    remove_over_nsamples: # done
        scatter: [vcf]
        run: ../tools/bcftools_isec_complement.cwl
        in:
            vcf: caller_vcfs
            exclude_vcf: nsamples/merged_nsamples_overN
            output_vcf_name: 
                valueFrom: "$(inputs.vcf.nameroot).nsamples.removed.vcf.gz"
        out:
            [complement_vcf]
        doc: "hard filter"
    pon_pileup: # done
        run: ../subworflows/get_base_counts_multi_sample.cwl
        in:
            reference: reference
            normal_bams: normal_bams
            merged_vcf: nsamples/merged_nsamples_underN
            mapq: mapq
            baseq: baseq
            pon_final_name: pon_final_name
        out:
            [total_counts]
    call_R_fisher: # need to do
        run: ../tools/normal_fisher.cwl
        scatter: [vcf]
        in:
            vcf: remove_over_nsamples/nsamples_filtered_vcf
            pon_total_counts: pon_pileup
            p-value: "1.293211e-09"
            caller: caller_prefix
        out:
            [pon_vcf, pon_filtered_vcf]
        doc: "hard and soft filter but only annotating hard filter with VEP"
    concat_vcfs: # need to do
        run: ../tools/concat_multisample_vcf.cwl
        in:
            vcfs: call_R_fisher/pon_vcf
            filtered_vcfs: call_R_fisher/pon_filtered_vcf
            final_prefix: $(caller_prefix).final.pon
        out:
            [pon_concatenated_vcf, pon_filtered_concatenated_vcf]
    