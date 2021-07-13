#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Using pyvcf count number of samples with variant in a multivcf. https://pyvcf.readthedocs.io/en/latest/FILTERS.html"
doc: "If the number of samples with variant is greater than (>) these are flagged and then removed"

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
            export varscan_vcf="$2"
            export vardict_vcf="$3"
            export pindel_vcf="$4"
            export i="$6"
            export b="$7"
            export c="$8"
            export T="$9"
            export oncoKB_curated="${10}"
            export p="${11}"
            export pan_myeloid="${12}"
            export blacklist="${13}"
            export segemental_duplications="${14}"
            export simple_repeats="${15}"
            export repeat_masker="${16}"


            LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate_CH_pd.R -v $mutect_vcf,$varscan_vcf,$vardict_vcf,$pindel_vcf \
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
    varscan_vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
    vardict_vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 3
    pindel_vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 4
    sample_name:
        type: string
        inputBinding:
            position: 5
    impact_annotation:
        type: File
        inputBinding:
            position: 6
    topmed_annotation:
        type: File
        inputBinding:
            position: 7
    cosmic_annotation:
        type: File
        inputBinding:
            position: 8
    tsg_annotation:
        type: File
        inputBinding:
            position: 9
    oncoKB_annotation:
        type: File
        inputBinding:
            position: 10
    pd_table_annotation:
        type: File
        inputBinding:
            position: 11
    panmyeloid_annotation:
        type: File
        inputBinding:
            position: 12
    blacklist_annotation:
        type: File
        inputBinding:
            position: 13
    segemental_duplications_annotation:
        type: File
        inputBinding:
            position: 14
    simple_repeats_annotation:
        type: File
        inputBinding:
            position: 15
    repeat_masker_annotation:
        type: File
        inputBinding:
            position: 16
    
outputs:
    final_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).final.tsv"
    column_check:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).columns.txt"
