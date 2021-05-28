#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'extract umis from bam'
baseCommand: ["/bin/bash", "extract_umis_helper.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'extract_umis_helper.sh'
        entry: |
            set -eo pipefail

            if [[ "$5" == 1 ]]; then
                /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA ZB --single-tag RX --input $1 --read-structure $2 $3 --output $4
            else
                /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA --single-tag RX --input $1 --read-structure $2 $3 --output $4
            fi

arguments:
    - output:
      valueFrom: $(runtime.outdir)/umi_extracted.bam
      position: 3
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    read_structure:
        type: string[]
        inputBinding:
            position: 2
    UMI_paired:
        type: int
        default: 1
        inputBinding:
            position: 5
outputs:
    umi_extracted_bam:
        type: File
        outputBinding:
            glob: "umi_extracted.bam"
