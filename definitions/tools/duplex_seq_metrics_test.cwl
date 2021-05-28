#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'collect duplex seq metrics'

baseCommand: ["/bin/bash", "collect_duplex_metrics_helper.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 10000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'collect_duplex_metrics_helper.sh'
        entry: |
            set -eo pipefail
            PAIRED=false
            while getopts "i:b:d:o:p" opt; do
                case "$opt" in
                    i)
                        INPUT="$OPTARG"
                        ;;
                    b)
                        INTERVALS="$OPTARG"
                        ;;
                    d)
                        DESCRIPTION="$OPTARG"
                        ;;
                    o)
                        OUTPUT="$OPTARG"
                        ;;
                    p)
                        PAIRED="$OPTARG"
                        ;;
                esac
            done

            if [[ $PAIRED == 1 ]]; then
                if [[ -z "$DESCRIPTION" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input "$INPUT" --description "$DESCRIPTION" --output "$OUTPUT"
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input "$INPUT" --intervals "$INTERVALS" --description "$DESCRIPTION" --output "$OUTPUT"
                fi
            else
                echo "Sample not UMI Paired" > "$OUTPUT"
            fi

arguments:
    ["-o", { valueFrom: "$(runtime.outdir)/duplex_seq.metrics"} ]

inputs:
    bam:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    intervals:
        type: File?
        inputBinding:
            prefix: "-b"
            position: 2
    description:
        type: string
        inputBinding:
            prefix: "-d"
            position: 3
    UMI_paired:
        type: int
        default: 1
        inputBinding:
            position: 5
            prefix: "-p"
outputs:
    duplex_seq_metrics:
        type: File[]
        outputBinding:
            glob: "duplex_seq.metrics.*"
