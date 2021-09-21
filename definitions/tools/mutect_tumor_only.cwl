#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Mutect2 (GATK 4)"

baseCommand: ["/bin/bash", "Mutect2.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'Mutect2.sh'
        entry: |
            set -o pipefail
            set -o errexit

            OUTPUT="$1"
            REF="$2"
            export tumor_bam="$3"
            BAM="$3"
            INTERVAL_LIST="$4"

            TUMOR=`samtools view -H $tumor_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
            echo -n $TUMOR > sampleName.txt

            /gatk/gatk Mutect2 --native-pair-hmm-threads 28 -R $REF -L $INTERVAL_LIST -I $BAM --max-reads-per-alignment-start 0 -O $OUTPUT
            /gatk/gatk FilterMutectCalls -V mutect.vcf.gz --reference $REF -O mutect.filtered.vcf.gz

arguments:
    - position: 1
      valueFrom: mutect.vcf.gz

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: 2
    tumor_bam:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.bai]
    interval_list:
        type: File
        inputBinding:
            position: 4

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "mutect.filtered.vcf.gz"
        secondaryFiles: [.tbi]
    tumor_sample_name:
        type: string
        outputBinding:
            glob: sampleName.txt
            loadContents: true
            outputEval: $(self[0].contents)
