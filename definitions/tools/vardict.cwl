#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "vardict.sh"]

requirements:
    - class: DockerRequirement
      dockerPull: "kboltonlab/vardictjava:1.0"
    - class: ResourceRequirement
      ramMin: 20000
      coresMin: 1
      tmpdirMin: 10000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "vardict.sh"
        entry: |
          #!/bin/bash

          set -eou pipefail
          
          REF="$1"
          AF_THR="$2"
          tumor_bam="$3"
          tumor_sample_name="$4"
          bed="$5"
          normal_bam="$6"
          normal_sample_name="$7"
          out="$8"

          /opt/VarDictJava/build/install/VarDict/bin/VarDict -G $REF -f $AF_THR -N $tumor_sample_name -b "$tumor_bam|$normal_bam" -c 1 -S 2 -E 3 -g 4 $bed | /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -N "$tumor_sample_name|$normal_sample_name" -f $AF_THR > $out

          bgzip $out && tabix $out.gz


inputs:
    reference:
        type: File
        inputBinding:
            position: 1
    af_threshold:
        type: float?
        inputBinding:
            position: 2
        default: 0.005
        doc: "minimum allele frequency"
    tumor_bam:
        type: File
        inputBinding:
            position: 3
    tumor_sample_name:
        type: string
        inputBinding:
            position: 4
        default: 'TUMOR'
    interval_list:
        type: File
        inputBinding:
            position: 5
        doc: "required bed file by VardictJava"
    normal_bam:
        type: File
        inputBinding:
            position: 6
    normal_sample_name:
        type: string
        inputBinding:
            position: 7
        default: 'NORMAL'
    output_name:
        type: string?
        inputBinding:
            position: 8
        default: 'vardict.vcf'

outputs:
    vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_name).gz
        secondaryFiles: [.tbi]

