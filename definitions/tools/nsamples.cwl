#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Using pyvcf count number of samples with variant in a multivcf. https://pyvcf.readthedocs.io/en/latest/FILTERS.html"
doc: "If the number of samples with variant is greater than (>) these are flagged and then removed"

baseCommand: ["/bin/bash", "combine.sh"]

arguments: [
    { position: 4, valueFrom: "$(inputs.nsamples+1)orGreaterRemoved" },
    { position: 5, valueFrom: "$(inputs.nsamples+1)orGreater" }
]

requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "kboltonlab/pyvcf:1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'combine.sh'
        entry: |
            set -eou pipefail
            set -o errexit

            ## the filter_nsamples.py python script has to be in the current directory
            echo "
            import vcf.filters

            class NSamples(vcf.filters.Base):
                'Filter sites by N samples'

                name = 'ns'

                @classmethod
                def customize_parser(self, parser):
                    parser.add_argument('--n-samples', type=int, default=5,
                            help='Filter out sites with greater n samples')

                def __init__(self, args):
                    self.threshold = args.n_samples

                def __call__(self, record):
                    if record.num_called > self.threshold:
                        return record
            " > filter_nsamples.py
            chmod u+x filter_nsamples.py

            export output="$1"
            export input="$2"
            export nsamples="$3"
            export keep="$4"
            export remove="$5"
            
            
            /opt/conda/bin/vcf_filter.py --local-script filter_nsamples.py --output "$output.vcf" "$input" ns --n-samples "$nsamples"

            # need to install bcftools, tabix, bgzip
            bgzip "$output.vcf" && tabix "$output.vcf.gz"
            filter_string="FILTER~\"ns$nsamples\""
            ## this is kept for MSK pileup
            bcftools filter -e $filter_string -Ov -o "$output.$keep.vcf" "$output.vcf.gz"
            ## this is kept to intersect-complement individual vcf files
            bcftools filter -i $filter_string -Oz -o "$output.$remove.vcf.gz" "$output.vcf.gz" && tabix "$output.$remove.vcf.gz"

inputs:
    output_name:
        type: string
        inputBinding:
            position: 1
    merged_vcf:
        type: File
        inputBinding:
            position: 2
    nsamples:
        type: int
        default: 5
        inputBinding:
            position: 3
    
outputs:
    merged_nsamples:
        type: File
        outputBinding:
            glob: "$(inputs.output_name).vcf.gz"
        secondaryFiles: [.tbi]
    merged_nsamples_underN:
        type: File
        outputBinding:
            glob: "$(inputs.output_name).$(inputs.nsamples+1)orGreaterRemoved.vcf"
        secondaryFiles: [.tbi]
    merged_nsamples_overN:
        type: File
        outputBinding:
            glob: "$(inputs.output_name).$(inputs.nsamples+1)orGreater.vcf.gz"
        secondaryFiles: [.tbi]

