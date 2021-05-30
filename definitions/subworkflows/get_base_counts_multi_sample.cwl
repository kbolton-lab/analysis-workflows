#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "MSK Pileup workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: InlineJavascriptRequirement
    - class: StepInputExpressionRequirement

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    normal_bams:
        type: File[]
        secondaryFiles: [.bai]
    vcf:
        type: File
        doc: "Must be in vcf format, not .gz.  Usually after filtered from Nsamples"
    mapq:
        type: int?
        default: 5
    baseq:
        type: int?
        default: 5
    pon_final_name:
        type: string?
        default: "pon.all"
        doc: "prefix final name for total pileup, usually 'pon.all' + caller name"
outputs:
    pileups:
        type: File[]
        outputSource: get_pileup_counts/pileup
    pileup_counts:
        type: File[]
        outputSource: get_pileup_counts/pileup_counts
    total_counts:
        type: File
        outputSource: get_total/total_counts
        secondaryFiles: [.tbi]
steps:
    get_pileup_counts:
        scatter: [bam]
        run: ../tools/msk_get_base_counts.cwl
        in:
            bam: normal_bams
            sample_name: 
                valueFrom: $(inputs.bam.nameroot)
            reference: reference
            vcf: vcf
            baseq: baseq
            mapq: mapq
        out:
            [pileup, pileup_counts]
    get_total:
        run: ../tools/total_msk_get_base_counts.cwl
        in:
            counts: get_pileup_counts/pileup_counts
            final_name: pon_final_name
        out:
            [total_counts]
    