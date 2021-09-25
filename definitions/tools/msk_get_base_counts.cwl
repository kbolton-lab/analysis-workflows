#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "msk_pileup.sh"]

arguments: [
    { position: 7, valueFrom: $(runtime.cores) }
]

requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 16
      ramMin: 128000
    - class: DockerRequirement
      dockerPull: "kboltonlab/msk_getbasecounts:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "msk_pileup.sh"
        entry: |

          set -eou pipefail

          export sample_name="$2"
          export bam_path="$3"
          # optionally can get samples name from but this will be different and maybe cannot outputBinding
          # export sample_name=`samtools view -H test.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`
          echo "$sample_name:$bam_path"

          /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta "$1" --bam "$sample_name:$bam_path"  --vcf "$4" --output "$sample_name.pileup.vcf" --maq "$5" --baq "$6" --thread 16;

          bgzip $sample_name.pileup.vcf && tabix $sample_name.pileup.vcf.gz

          bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%RD]\t[%AD]\n' $sample_name.pileup.vcf.gz > $sample_name.pileup.txt
inputs:
    reference:
        type:
            - string
            - File
        inputBinding:
            position: 1
    bam:
        type: string
        inputBinding:
            position: 3
        doc: "Normal bam path to storage1"
    sample_name:
        type: string
        inputBinding:
            position: 2
    vcf:
        type: File
        inputBinding:
            position: 4
        doc: "VCF file to perform counts"
    mapq:
        type: int
        inputBinding:
            position: 5
    baseq:
        type: int
        inputBinding:
            position: 6

outputs:
    pileup:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).pileup.vcf.gz"
    pileup_counts:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).pileup.txt"
