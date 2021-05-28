#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Total MSK"

baseCommand: ["/bin/bash", "combine.sh"]

requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 2000
    - class: DockerRequirement
      dockerPull: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'combine.sh'
        entry: |

            set -eou pipefail

            export counts="$1"
            export output="$2"
            IFS=' ' read -r array <<< "$counts"

            awk '{chrom[FNR]=$1; pos[FNR]=$2; ref[FNR]=$3; alt[FNR]=$4; b[FNR]+=$5; c[FNR]+=$6}END{print "Chrom\tPOS\tREF\tALT\tTotal_Ref_Depth\tTotal_Alt_Counts\tTotal_Alt_Freq"; for(i=1;i<=FNR;i++) {
                if ((c[i]+b[i])==0) printf "%s",chrom[i]"\t"pos[i]"\t"ref[i]"\t"alt[i]"\t"b[i]"\t"c[i]"\t0.0\n";
                else printf "%s",chrom[i]"\t"pos[i]"\t"ref[i]"\t"alt[i]"\t"b[i]"\t"c[i]"\t"(c[i]/(c[i]+b[i]))"\n";}
            }' `printf '%s\n' $array` > "$output";
            
            tail -n +2 "$output" | bgzip -c > "$output.gz" && tabix -f -s1 -b2 -e2 "$output.gz";

inputs:
    counts:
        type: 
            type: array
            items: File
        inputBinding:
            position: 1
            itemSeparator: " " 
    final_name:
        type: string
        default: "all.counts"
        inputBinding:
            position: 2

outputs:
    total_counts:
        type: File
        outputBinding:
            glob: "$(inputs.final_name).gz"
        secondaryFiles: [.tbi]