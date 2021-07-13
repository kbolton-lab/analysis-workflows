#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "gnomad parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
inputs:
    gnomad_AF_only:
        type: File
        secondaryFiles: [.tbi]
    gnomad_AF:
        type: float
        default: .005
    info_af_filter_string: 
        type: string
        default: "INFO/AF>"
    filter_flag:
        type:
            type: enum
            symbols: ["include", "exclude"]
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    # output_prefix: 
    #     type: string
    #     default: "caller"
outputs:
    gnomad_exclude:
        type: File
        outputSource: filter_gnomAD_exclude/filtered_vcf
        secondaryFiles: [.tbi]
    normalized_gnomad_exclude:
        type: File
        outputSource: index_exclude_norm/indexed_vcf
        secondaryFiles: [.tbi]
    # filtered_vcf:
    #     type: File
    #     outputSource: isec_complement_gnomAD/complement_vcf
    

steps:
    filter_gnomAD_exclude:
        run: ../tools/bcftools_filter.cwl
        in:
            vcf: gnomad_AF_only
            filter_string: 
                source: ["#info_af_filter_string", "#gnomad_AF"]
                valueFrom: ${ return self[0] + self[1].toString(); }
            filter_flag: filter_flag
            output_type:
                default: "z"
            output_vcf_name: 
                valueFrom: "gnomad.AF.exclude.vcf.gz"
        out: [filtered_vcf]
        doc: "this filter's the gnomAD_af_only file based on gnomAD POPAF threshold, it is what should be excluded if our calls have it since above threshold"
    normalize_gnomAD_exclude:
        run: ../tools/bcftools_norm.cwl
        in:
            vcf: filter_gnomAD_exclude/filtered_vcf
            output_type:
                default: "z"
            output_vcf_name: 
                valueFrom: "gnomad.AF.exclude.norm.vcf.gz"
        out: [normalized_vcf]
        doc: "for isec to work on multiallelics, the file must be normalized"
    index_exclude_norm:
        run: ../tools/index_vcf.cwl
        in:
            vcf: normalize_gnomAD_exclude/normalized_vcf
        out:
            [indexed_vcf]
    # isec_complement_gnomAD:
    #     run: ../tools/bcftools_isec_complement.cwl
    #     in:
    #         vcf: vcf
    #         exclude_vcf: index_exclude_norm/indexed_vcf
    #         output_vcf_name: 
    #             source: "#output_prefix"
    #             valueFrom: "$(self).gnomAD_AF_filter.vcf"
    #         output_type: 
    #             default: "v"
    #     doc: "both input files should be bgzipped and indexed"
    #     out:
    #         [complement_vcf]
