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
            type: enum
            symbols: ["include", "exclude"]
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    nsamples:
        type: int
        default: 5
        doc: |
            The threshold for number of samples a mutation can have to filter out artefacts and sequencing errors
            Can be found with qbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE) in R. If number of samples with variant > nsamples these are removed
    mapq:
        type: int?
        default: 0
    baseq:
        type: int?
        default: 0
outputs:
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
    remove_over_nsamples:
        type: File[]
        outputSource: remove_over_nsamples/complement_vcf

steps:
    merge:
        run: ../tools/bcftools_merge.cwl
        in:
            vcfs: caller_vcfs
            merge_method: merge_method
            output_vcf_name: 
                source: "#caller_prefix"
                valueFrom: "$(self).merged.vcf.gz"
            output_type: 
                default: "z"
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
            filter_flag: filter_flag
            output_prefix: caller_prefix
        out:
            [filtered_vcf]
    gnomAD_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: gnomAD_filter/filtered_vcf
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
    remove_over_nsamples: # need to do
        scatter: [vcf]
        run: ../tools/bcftools_isec_complement.cwl
        in:
            vcf: caller_vcfs
            exclude_vcf: nsamples/merged_nsamples_overN
            output_vcf_name: 
                valueFrom: "$(inputs.vcf.nameroot).nsamples.removed.vcf.gz"
        out:
            [complement_vcf]
    