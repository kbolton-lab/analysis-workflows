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

            export tumor_bam="$3"
            export normal_bam="$5"

            TUMOR=`samtools view -H $tumor_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
            echo -n $TUMOR > sampleName.txt
            echo -n $normal_bam >> sampleName.txt
            # echo "{'tumor_sample_name': $TUMOR}"
            if [ -z "${normal_bam}" ]; then
                /gatk/gatk Mutect2 --java-options "-Xmx20g" --native-pair-hmm-threads 28 -R $2 -L $4 -I $3 --max-reads-per-alignment-start 0 -O $1
                /gatk/gatk FilterMutectCalls -V mutect.vcf.gz --reference $2 -O mutect.filtered.vcf.gz
            else
                NORMAL=`samtools view -H $normal_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
                /gatk/gatk Mutect2 --java-options "-Xmx20g" --native-pair-hmm-threads 28 -O $1 -R $2 -I $3 -tumor "$TUMOR" -I $5 -normal "$NORMAL" -L $4 --max-reads-per-alignment-start 0 #Running Mutect2.
                /gatk/gatk FilterMutectCalls -R $2 -V mutect.vcf.gz -O mutect.filtered.vcf.gz #Running FilterMutectCalls on the output vcf.
            fi

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
    normal_bam:
        type: File?
        inputBinding:
            position: 5
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
    # tumor_sample_name: string
    tumor_sample_name:
        type: string
        outputBinding:
            glob: sampleName.txt
            loadContents: true
            outputEval: $(self[0].contents)
