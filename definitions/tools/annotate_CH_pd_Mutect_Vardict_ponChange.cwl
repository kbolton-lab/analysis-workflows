#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "R annotation and filtering"
doc: "Run combining of the variants and annotate_PD"

baseCommand: ["/bin/bash", "annotate.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 16000
      tmpdirMin: 100000
      coresMin: 6
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "kboltonlab/annotate_wes_ch:3.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'annotate.sh'
        entry: |
            
            set -eou pipefail
        
            export mutect_vcf="$1";
            shift;
            export vardict_vcf="$1";
            shift;
            export sample="$1";
            shift;
            i="$1"
            shift;
            b="$1"
            shift;
            c="$1"
            shift;
            T="$1"
            shift;
            oncoKB_curated="$1"
            p="$2"
            pan_myeloid="$3"
            blacklist="$4"
            segemental_duplications="$5"
            simple_repeats="$6"
            repeat_masker="$7"
            intervals="$8"

            target_length=`grep -v "^@" $intervals | awk '{print $1,$2,$3}' | awk '$4=$3-$2 {sum+=$4}END {print sum}'`

            LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate_CH_pd_docker3_Mutect_Vardict_ponChange.R -v $mutect_vcf,$vardict_vcf \
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
                                                                       --repeat-masker $repeat_masker \
                                                                       --target-length $target_length


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
    intervals:
        type: File
        inputBinding:
            position: 15
    
outputs:
    final_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).final.tsv"
    column_check:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).columns.txt"
