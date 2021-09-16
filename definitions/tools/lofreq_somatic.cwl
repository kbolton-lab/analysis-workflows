#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "lofreq v2.1.3.1 call"
baseCommand: ["/bin/bash", "lofreq_helper.sh"]
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

            if [ $# -lt 3 ]
            then
                echo "Usage: $0 [TUMOR_BAM] REFERENCE] [roi_bed?]"
                exit 1
            fi

            TUMOR_BAM="$1"
            NORMAL_BAM="$2"
            REFERENCE="$3"
            OUTPUT="$4"

            if [ -z ${3+x} ]; then
                #run without ROI
                /opt/lofreq/bin/lofreq somatic -n $NORMAL_BAM -t $TUMOR_BAM -f $REFERENCE -o $OUTPUT
            else
                ROI_BED="$4"
                /opt/lofreq/bin/lofreq somatic -n $NORMAL_BAM -t $TUMOR_BAM -f $REFERENCE -l $ROI_BED -o $OUTPUT
            fi
            tabix $OUTPUTsomatic_final.snvs.vcf.gz
            tabix $OUTPUTsomatic_final.indels.vcf.gz
            bcftools concat -a $OUTPUTsomatic_final.snvs.vcf.gz $OUTPUTsomatic_final.indels.vcf.gz > $OUTPUT.vcf
            bgzip $OUTPUT.vcf && tabix $OUTPUT.vcf.gz

inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.bai]
    reference:
        type: File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: 3
    interval_list:
        type: File?
        inputBinding:
            position: 5
    output_name:
        type: string?
        inputBinding:
            position: 4
        default: "lofreq"

outputs:
    vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_name).vcf.gz
        secondaryFiles: [.tbi]
