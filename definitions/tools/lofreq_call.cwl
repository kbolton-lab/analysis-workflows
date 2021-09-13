#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "lofreq v2.1.3.1 call"
baseCommand: "lofreq_helper.sh"
requirements:
    - class: ResourceRequirement
      ramMin: 12000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: "kboltonlab/lofreq:latest"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'lofreq_helper.sh'
        entry: |
            #!/bin/bash

            set -o errexit
            set -o nounset

            if [ $# -lt 2 ]
            then
                echo "Usage: $0 [TUMOR_BAM] REFERENCE] [roi_bed?]"
                exit 1
            fi

            TUMOR_BAM="$1"
            REFERENCE="$2"
            OUTPUT="$3"

            if [ -z ${3+x} ]; then
                #run without ROI
                lofreq call -A -B -f $REFERENCE --call-indels -o $OUTPUT $TUMOR_BAM
            else
                ROI_BED="$3"
                lofreq call -A -B -f $REFERENCE --call-indels --bed $ROI_BED -o $OUTPUT $TUMOR_BAM
            fi
            bgzip $OUTPUT && tabix $OUTPUT.gz
inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: 2
    roi_bed:
        type: File?
        inputBinding:
            position: 3
    output_name:
        type: string?
        inputBinding:
            position: 4
        default: "lofreq.vcf"

outputs:
    vcfs:
        type: File
        outputBinding:
            glob: $(inputs.output_name).gz
        secondaryFiles: [.tbi]