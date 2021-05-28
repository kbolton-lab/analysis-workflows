#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/gatk-cwl:3.6.0
arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "Mutect2,Vardict,Varscan2,Pindel",
     "-o", { valueFrom: "$(runtime.outdir)/" + "$(inputs.tumor_sample_name)" + ".combined.vcf.gz" }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    input_vcfs:
        type: File[]
        inputBinding:
        prefix: -C=
        itemSeparator: " "
        separate: false
        position: 2
    mutect_vcf:
        type: File
        inputBinding:
            prefix: "--variant:Mutect2"
            position: 2
        secondaryFiles: [.tbi]
    varscan_vcf:
        type: File
        inputBinding:
            prefix: "--variant:Varscan2"
            position: 3
        secondaryFiles: [.tbi]
    vardict:
        type: File
        inputBinding:
            prefix: "--variant:Vvardict"
            position: 4
        secondaryFiles: [.tbi]
    pindel_vcf:
        type: File
        inputBinding:
            prefix: "--variant:Pindel"
            position: 5
        secondaryFiles: [.tbi]
    tumor_sample_name :
        type: string
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: 
                valueFrom: "$(inputs.tumor_sample_name)" + ".combined.vcf.gz"

