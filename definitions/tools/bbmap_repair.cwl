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
arguments: [ "in1=", { valueFrom: $(inputs.fastq1) }, "in2=", { valueFrom: $(inputs.fastq2) }]
inputs:
    fastq1:
        type: File
    fastq2:
        type: File
outputs:
    fastqs:
        type: File[]
        outputBinding:
            glob: "*.fixed.fastq.gz"
    fastq1:
        type: File
        outputBinding:
            glob: "R1.fixed.fastq.gz"
    fastq2:
        type: File
        outputBinding:
            glob: "R2.fixed.fastq.gz"
