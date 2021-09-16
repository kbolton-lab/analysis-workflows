#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "AWK: Keep UMI of Certain Length"
baseCommand: ["/bin/bash", "filter_umi_length_helper.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
    - class: SchemaDefRequirement
      types:
        - $import: ../types/sequence_data.yml
    - class: InitialWorkDirRequirement
      listing:
        - entryname: "filter_umi_length_helper.sh"
          entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            FASTQ_ONE="$1"
            FASTQ_TWO="$2"
            UMI_LENGTH="$3"
            NAME="$4"

            zcat $FASTQ_ONE | awk -v regex="AACCGCCAGGAGT" length="$UMI_LENGTH" 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; split(seq,a,regex); if (length(a[1]) == length) {print header, seq, qheader, qseq}}' > $NAME.fastq.gz
arguments:
    - {valueFrom: "$(inputs.sequence.sequence.fastq1)"}
    - {valueFrom: "$(inputs.sequence.sequence.fastq2)"}
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    umi_length:
        type: int
        inputBinding:
            position: 3
    output_fastq_name:
        type: string
        default: "R1_filtered"
outputs:
    fastq1:
        type: File
        outputBinding:
            glob: "R1_filtered.fastq.gz"
    fastq2:
        type: File
        outputBinding:
            glob: $(inputs.sequence.sequence.fastq2)
