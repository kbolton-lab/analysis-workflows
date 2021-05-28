#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'group reads by umi'
baseCommand: ["/bin/bash", "group_reads_helper.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'group_reads_helper.sh'
        entry: |
            set -eo pipefail

            if [[ "$3" == 1 ]]; then
                /usr/local/bin/fgbio GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $1 --output $2
            else
                /usr/local/bin/fgbio GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $1 --output $2
            fi

arguments:
    - output:
      valueFrom: $(runtime.outdir)/umi_grouped.bam
      position: 2

inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    UMI_paired:
        type: int
        default: 1
        inputBinding:
            position: 3
outputs:
    grouped_bam:
        type: File
        outputBinding:
            glob: "umi_grouped.bam"
