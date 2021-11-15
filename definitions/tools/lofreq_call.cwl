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

            if [ $# -lt 4 ]
            then
                echo "Usage: $0 [TUMOR_BAM] REFERENCE] [OUTPUT_NAME] [MIN_VAF] [roi_bed?]"
                exit 1
            fi

            TUMOR_BAM="$1"
            REFERENCE="$2"
            OUTPUT="$3"
            MIN_VAF="$4"

            if [ -z ${4+x} ]; then
                #run without ROI
                /opt/lofreq/bin/lofreq indelqual --dindel -f $REFERENCE -o output.indel.bam $TUMOR_BAM
                /opt/lofreq/bin/lofreq call -A -B -f $REFERENCE --call-indels -o lofreq_pass.vcf output.indel.bam --force-overwrite

                /opt/lofreq/bin/lofreq call --no-default-filter -B -a 1 -b 1 -f $REFERENCE --call-indels -o lofreq_call.vcf output.indel.bam --force-overwrite
                /opt/lofreq/bin/lofreq filter -i lofreq_call.vcf -o lofreq_call.filtered.vcf -v 5 -a ${MIN_VAF} -A 0.9 --sb-incl-indels --print-all
            else
                ROI_BED="$5"
                /opt/lofreq/bin/lofreq indelqual --dindel -f $REFERENCE -o output.indel.bam $TUMOR_BAM
                /opt/lofreq/bin/lofreq call -A -B -f $REFERENCE --call-indels --bed $ROI_BED -o lofreq_pass.vcf output.indel.bam --force-overwrite

                /opt/lofreq/bin/lofreq call --no-default-filter -B -a 1 -b 1 -l $ROI_BED -f $REFERENCE --call-indels -o lofreq_call.vcf output.indel.bam --force-overwrite
                /opt/lofreq/bin/lofreq filter -i lofreq_call.vcf -o lofreq_call.filtered.vcf -v 5 -a ${MIN_VAF} -A 0.9 --sb-incl-indels --print-all
            fi

            printf "##FILTER=<ID=CALL,Description=\"A variant that was called by Lofreq's Caller without any filters\">" > lofreq.header;
            cat lofreq_call.filtered.vcf | sed 's/PASS/CALL/g' > call_to_pass.vcf
            bgzip call_to_pass.vcf && tabix call_to_pass.vcf.gz
            bgzip lofreq_pass.vcf && tabix lofreq_pass.vcf.gz
            bcftools annotate --threads 32 -a lofreq_pass.vcf.gz -h lofreq.header -c FILTER call_to_pass.vcf.gz -Oz -o $OUTPUT.gz

            tabix $OUTPUT.gz
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
    interval_list:
        type: File?
        inputBinding:
            position: 5
    min_var_freq:
        type: float?
        inputBinding:
            position: 4
        default: 0.005
    output_name:
        type: string?
        inputBinding:
            position: 3
        default: "lofreq.vcf"

outputs:
    vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_name).gz
        secondaryFiles: [.tbi]
