#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "BBMap Repair Fastq1 with Fastq2"
baseCommand: ["repair.sh", "repair=t", "overwrite=true", "interleaved=false", "outs=singletons.fq", "out1=R1.fixed.fastq.gz", "out2=R2.fixed.fastq.gz"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
arguments: [
    - {valueFrom: "in1=$(inputs.fastq1)"}
    - {valueFrom: "in2=$(inputs.fastq2)"}
]
inputs:
    fastq1:
        type: File[]
    fastq2:
        type: File[]
outputs:
    fastq1:
        type: File
        outputBinding:
            glob: "R1.fixed.fastq.gz"
    fastq2:
        type: File
        outputBinding:
            glob: "R2.fixed.fastq.gz"
