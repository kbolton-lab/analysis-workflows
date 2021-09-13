#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "R annotation and filtering"
doc: "Run combining of the variants and annotate_PD"

baseCommand: ["/bin/bash", "annotate.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
      coresMin: 12
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "kboltonlab/annotate_wes_ch:3.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'annotate.sh'
        entry: |
            
            set -eou pipefail
        
            export mutect_vcf="$1"
            export vardict_vcf="$2"
            export i="$4"
            export b="$5"
            export c="$6"
            export T="$7"
            export oncoKB_curated="${8}"
            export p="${9}"
            export pan_myeloid="${10}"
            export blacklist="${11}"
            export segemental_duplications="${12}"
            export simple_repeats="${13}"
            export repeat_masker="${14}"


            LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate_CH_pd_docker3_Mutect_Vardict.R -v $mutect_vcf,$vardict_vcf \
                                                                       -i $i \
                                                                       -b $b \
                                                                       -c $c \
                                                                       -T $T \
                                                                       --oncoKB-curated $oncoKB_curated \
                                                                       -p $p \
                                                                       --pan-myeloid $pan_myeloid \
                                                                       --blacklist $blacklist \
                                                                       --segemental-duplications $segemental_duplications \
                                                                       --simple-repeats $simple_repeats \
                                                                       --repeat-masker $repeat_masker


inputs:
    mutect_vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 1
    vardict_vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
    sample_name:
        type: string
        inputBinding:
            position: 3
    impact_annotation:
        type: File
        inputBinding:
            position: 4
    topmed_annotation:
        type: File
        inputBinding:
            position: 5
    cosmic_annotation:
        type: File
        inputBinding:
            position: 6
    tsg_annotation:
        type: File
        inputBinding:
            position: 7
    oncoKB_annotation:
        type: File
        inputBinding:
            position: 8
    pd_table_annotation:
        type: File
        inputBinding:
            position: 9
    panmyeloid_annotation:
        type: File
        inputBinding:
            position: 10
    blacklist_annotation:
        type: File
        inputBinding:
            position: 11
    segemental_duplications_annotation:
        type: File
        inputBinding:
            position: 12
    simple_repeats_annotation:
        type: File
        inputBinding:
            position: 13
    repeat_masker_annotation:
        type: File
        inputBinding:
            position: 14
    
outputs:
    final_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).final.tsv"
    column_check:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).columns.txt"
