#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/bash', 'cat_all_helper.sh']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
    - class: ResourceRequirement
      ramMin: 16000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'cat_all_helper.sh'
        entry: |
            set -eou pipefail

            /bin/cat "$@" > cat.temp;
            /bin/grep "ChrID" cat.temp > /dev/stdout
stdout: "all_region_pindel.head"
inputs:
    region_pindel_outs:
        type: File[]
        inputBinding:
            position: -1
outputs:
    all_region_pindel_head:
        type: stdout
