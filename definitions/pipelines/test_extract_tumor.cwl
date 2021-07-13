#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
          - $import: ../types/vep_custom_annotation.yml
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
label: "CH_exome: exome alignment and somatic variant detection"
doc: |
  CH_exome_cram is designed to perform processing of mutant/wildtype H.sapiens
  exome sequencing data (in CRAM format from UKBB) for low VAF variants. It features 4 caller variant
  detection, Nsamples and PoN filter, and vep style annotations.

  example input file = analysis_workflows/example_data/CH_exome_cram.yaml

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
        label: "reference: Reference fasta file for a desired assembly"
        doc: |
          reference contains the nucleotide sequence for a given assembly (hg37, hg38, etc.)
          in fasta format for the entire genome. This is what reads will be aligned to.
          Appropriate files can be found on ensembl at https://ensembl.org/info/data/ftp/index.html
          When providing the reference secondary files corresponding to reference indices must be
          located in the same directory as the reference itself. These files can be created with
          samtools index, bwa index, and picard CreateSequenceDictionary.
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai, .bai]
    tumor_name:
        type: string?
        default: 'tumor'
        label: "tumor_name: String specifying the name of the MT sample"
        doc: |
          tumor_name provides a string for what the MT sample will be referred to in the various
          outputs, for example the VCF files.  This should be the UKBB eid not the sample id???
    normal_name:
        type: string?
        default: 'normal'
        label: "normal_name: String specifying the name of the WT unmatched sample"
        doc: |
          normal_name provides a string for what the WT sample will be referred to in the various
          outputs, for example the VCF files.  This should be the UKBB eid not the sample id???
    tumor_sample_name:
        type: string
        doc: |
          this is the sample name, i.e. for UKBB it is UKB_1626332_232918706.  Should be given in 
          process_inputs.pl by: if($id->sample eq $tumor_sample) { $tumor_sample_name = $sm
    normal_sample_name:
        type: string
        doc: |
          this is the sample name, i.e. for UKBB it is UKB_1626332_232918706.  Should be given in 
          process_inputs.pl by: elsif($id->sample eq $normal_sample) { $normal_sample_name = $sm
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the capture reagent.
          Bed files with this information can be converted to interval_lists with Picard BedToIntervalList.
          In general for a WES exome reagent bait_intervals and target_intervals are the same. UKBB is xgen_plus_spikein.GRCh38.intervals
    target_interval_padding:
        type: int
        label: "target_interval_padding: number of bp flanking each target region in which to allow variant calls"
        doc: |
          The effective coverage of capture products generally extends out beyond the actual regions
          targeted. This parameter allows variants to be called in these wingspan regions, extending
          this many base pairs from each side of the target regions.
        default: 100
    qc_minimum_mapping_quality:
        type: int?
        default: 0
    qc_minimum_base_quality:
        type: int?
        default: 0
    # nsamples:
    #     type: int?
    #     default: 15
    #     doc: |
    #         The minimum alloted threshold for number of samples a mutation can have to filter out artefacts and sequencing errors
    #         Can be found with qbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE) in R. 15 for UKBB Batch qbinom(0.95, 1000, .01)
    scatter_count:
        type: int
        default: 2
        doc: "scatters each supported variant detector (varscan, pindel, mutect) into this many parallel jobs"
    mutect_artifact_detection_mode:
        type: boolean
        default: false
    mutect_max_alt_allele_in_normal_fraction:
        type: float?
    mutect_max_alt_alleles_in_normal_count:
        type: int?
    varscan_strand_filter:
        type: int?
        default: 0
    varscan_min_coverage:
        type: int?
        default: 4
    varscan_min_var_freq:
        type: float?
        default: 0.005
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_max_normal_freq:
        type: float?
    pindel_insert_size:
        type: int
        default: 400

 


outputs:
    tumor_bam:
        type: File
        outputSource: index_bam/indexed_bam
    mutect_unfiltered:
        type: File
        outputSource: mutect/unfiltered_vcf
    mutect_extract_tumor_index_file:
        type: File
        outputSource: mutect_extract_tumor_index/indexed_vcf
 

steps:
    index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: tumor_cram
        out:
            [indexed_cram]
    tumor_cram_to_bam:
        run: ../subworkflows/cram_to_bam_and_index.cwl
        in:
            reference: reference
            cram: index_cram/indexed_cram
        out:
            [bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: tumor_cram_to_bam/bam
        out:
            [indexed_bam] 
    mutect:
        run: ../subworkflows/mutect_normalize.cwl
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            normal_bam: normal_bam
            # interval_list: pad_target_intervals/expanded_interval_list
            interval_list: target_intervals
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    mutect_extract_tumor:
        run: ../tools/bcftools_extract_tumor.cwl
        in:
            vcf: mutect/unfiltered_vcf
            output_type:
                default: "z"
            tumor_sample_name: tumor_sample_name
            output_vcf_name: 
                valueFrom: |
                    $(inputs.tumor_sample_name).vcf.gz
        out:
            [tumor_only_vcf]
    mutect_extract_tumor_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: mutect_extract_tumor/tumor_only_vcf
        out:
            [indexed_vcf]