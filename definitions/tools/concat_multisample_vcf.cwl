#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Concat multisample for read in by vcfR"

baseCommand: ["/bin/bash", "concat.sh"]
arguments:
  - position: 4
    valueFrom: $(inputs.vcfs[0])
  - position: 5
    valueFrom: $(inputs.filtered_vcfs[0])


requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 2000
    - class: DockerRequirement
      dockerPull: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
    - class: InitialWorkDirRequirement 
      listing:
      - entryname: 'concat.sh'
        entry: |

            set -eou pipefail

            export vcfs="$1"
            export filtered_vcfs="$2"
            export prefix="$3"
            IFS=' ' read -r vcfs_array <<< "$vcfs"
            IFS=' ' read -r filtered_vcfs_array <<< "$filtered_vcfs"


            bcftools view -h "$4" > blanket.header
            sed -i 's/FORMAT\t.*/FORMAT\tFORMAT_STATS/' blanket.header

            awk '
            FNR==1 && NR!=1 { while (/^#/) getline; }
                1 {print}
            ' `printf '%s\n' $vcfs_array` > "$prefix.vcf"
            (cat blanket.header; bcftools view -H "$prefix.vcf") > "$prefix.vcf.tmp" && mv "$prefix.vcf.tmp" "$prefix.vcf"
            bcftools sort "$prefix.vcf" | bgzip -c > "$prefix.vcf.gz" && tabix "$prefix.vcf.gz"

            bcftools view -h "$5" > blanket.header
            sed -i 's/FORMAT\t.*/FORMAT\tFORMAT_STATS/' blanket.header
            awk '
            FNR==1 && NR!=1 { while (/^#/) getline; }
                1 {print}
            ' `printf '%s\n' $filtered_vcfs_array` > "$prefix.filtered.vcf"
            (cat blanket.header; bcftools view -H "$prefix.filtered.vcf") > "$prefix.filtered.vcf.tmp" && mv "$prefix.filtered.vcf.tmp" "$prefix.filtered.vcf"
            bcftools sort "$prefix.filtered.vcf" | bgzip -c > "$prefix.filtered.vcf.gz" && tabix "$prefix.filtered.vcf.gz"

inputs:
    vcfs:
        type: 
            type: array
            items: File
        inputBinding:
            position: 1
            itemSeparator: " " 
    filtered_vcfs:
        type: 
            type: array
            items: File
        inputBinding:
            position: 2
            itemSeparator: " "
    final_prefix:
        type: string?
        default: "caller.final.pon"
        inputBinding:
            position: 3

outputs:
    pon_concatenated_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.final_prefix).vcf.gz"
        secondaryFiles: [.tbi]
    pon_filtered_concatenated_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.final_prefix).filtered.vcf.gz"
        secondaryFiles: [.tbi]