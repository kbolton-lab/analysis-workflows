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
    tumor_crams:
        type: 
            type: array
            items: File
        secondaryFiles: [^.crai]
    # tumor_sample_names:
    #     type: string[]
    #     valueFrom: |
    #         ${
    #         var names = [];
    #         for( var i = 0; i < inputs.tumor_crams.length; i++) {
    #             names.push(inputs.tumor_crams[i].nameroot)
    #         }
    #         return names;
    #         }
    #     doc: "if this doesn't work then just get an array of names"
    # tumor_name:
    #     type: string?
    #     #default: $(inputs.tumor_crams.nameroot)
    #     default: 'tumor'
    #     label: "tumor_name: String specifying the name of the MT sample"
    #     doc: |
    #       tumor_name provides a string for what the MT sample will be referred to in the various
    #       outputs, for example the VCF files.
    normal_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    normal_sample_name:
        type: string?
        default: 'normal'
        label: "normal_sample_name: String specifying the name of the WT unmatched sample"
        doc: |
          normal_sample_name provides a string for what the WT sample will be referred to in the various
          outputs, for example the VCF files.
    pon_normal_bams:
        type: 
            type: array
            items: File
        secondaryFiles: [.bai,^.bai]
    #trimming:
    #    type:
    #        - ../types/trimming_options.yml#trimming_options
    #        - "null"
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        label: "bqsr_known_sites: One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
        doc: |
          Known polymorphic indels recommended by GATK for a variety of
          tools including the BaseRecalibrator. This is part of the GATK resource
          bundle available at http://www.broadinstitute.org/gatk/guide/article?id=1213
          File should be in vcf format, and tabix indexed.
    bqsr_intervals:
        type: string[]
        label: "bqsr_intervals: Array of strings specifying regions for base quality score recalibration"
        doc: |
          bqsr_intervals provides an array of genomic intervals for which to apply
          GATK base quality score recalibrations. Typically intervals are given
          for the entire chromosome (chr1, chr2, etc.), these names should match
          the format in the reference file.
    bait_intervals:
        type: File
        label: "bait_intervals: interval_list file of baits used in the sequencing experiment"
        doc: |
          bait_intervals is an interval_list corresponding to the baits used in sequencing reagent.
          These are essentially coordinates for regions you were able to design probes for in the reagent.
          Typically the reagent provider has this information available in bed format and it can be
          converted to an interval_list with Picard BedToIntervalList. Astrazeneca also maintains a repo
          of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the capture reagent.
          Bed files with this information can be converted to interval_lists with Picard BedToIntervalList.
          In general for a WES exome reagent bait_intervals and target_intervals are the same.
    target_interval_padding:
        type: int
        label: "target_interval_padding: number of bp flanking each target region in which to allow variant calls"
        doc: |
            The effective coverage of capture products generally extends out beyond the actual regions
            targeted. This parameter allows variants to be called in these wingspan regions, extending
            this many base pairs from each side of the target regions.
        default: 100
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        label: "per_base_intervals: additional intervals over which to summarize coverage/QC at a per-base resolution"
        doc: "per_base_intervals is a list of regions (in interval_list format) over which to summarize coverage/QC at a per-base resolution."
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        label: "per_target_intervals: additional intervals over which to summarize coverage/QC at a per-target resolution"
        doc: "per_target_intervals list of regions (in interval_list format) over which to summarize coverage/QC at a per-target resolution."
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
    qc_minimum_mapping_quality:
        type: int?
        default: 0
    qc_minimum_base_quality:
        type: int?
        default: 0
    nsamples:
        type: int?
        default: 15
        doc: |
            The minimum alloted threshold for number of samples a mutation can have to filter out artefacts and sequencing errors
            Can be found with qbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE) in R. 15 for UKBB Batch qbinom(0.95, 1000, .01)
    scatter_count:
        type: int
        default: 10
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
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.005
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_max_normal_freq:
        type: float?
    vardict_bed_targets:
        type: File
        doc: "VardictJava requires a bed file to run the variant calling, i.e. '-c 1 -S 2 -E 3 -g 4 /path/to/my.bed'"
    pindel_insert_size:
        type: int
        default: 400
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
        doc: "The set of alleles that gatk haplotype caller will use to force-call regardless of evidence"
    filter_docm_variants:
        type: boolean?
        default: true
    filter_somatic_llr_threshold:
        type: float
        default: 5
        doc: "Sets the stringency (log-likelihood ratio) used to filter out non-somatic variants.  Typical values are 10=high stringency, 5=normal, 3=low stringency. Low stringency may be desirable when read depths are low (as in WGS) or when tumor samples are impure."
    filter_somatic_llr_tumor_purity:
        type: float
        default: 1
        doc: "Sets the purity of the tumor used in the somatic llr filter, used to remove non-somatic variants. Probably only needs to be adjusted for low-purity (< 50%).  Range is 0 to 1"
    filter_somatic_llr_normal_contamination_rate:
        type: float
        default: 0
        doc: "Sets the fraction of tumor present in the normal sample (range 0 to 1), used in the somatic llr filter. Useful for heavily contaminated adjacent normals. Range is 0 to 1"
    vep_cache_dir:
        type:
            - string
            - Directory
        doc: "path to the vep cache directory, available at: https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#pre"
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    synonyms_file:
        type: File?
        doc: "synonyms_file allows the use of different chromosome identifiers in vep inputs or annotation files (cache, database, GFF, custom file, fasta). File should be tab-delimited with the primary identifier in column 1 and the synonym in column 2."
    annotate_coding_only:
        type: boolean?
        default: true
        doc: "if set to true, vep only returns consequences that fall in the coding regions of transcripts"
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
        doc: "configures how vep will annotate genomic features that each variant overlaps; for a detailed description of each option see https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_allele_gene_eg"
    filter_gnomADe_maximum_population_allele_frequency:
        type: float

    cle_vcf_filter:
        type: boolean
        default: false
    variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
        doc: "The names of one or more standard VCF fields or INFO fields to include in the output table"
    variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
        doc: "The name of a genotype field to include in the output table"
    vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
        doc: "VEP fields in final output"
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    
    normal_sample_name:
        type: string
        doc: "for varscan, pindel "
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "An optional VCF with variants that will be flagged as 'VALIDATED' if found in this pipeline's main output VCF"
    bed_targets:
        type: File
outputs:
    tumor_bam:
        type: File[]
        outputSource: tumor_cram_to_bam/bam
    tumor_mark_duplicates_metrics:
        type: File[]
        outputSource: tumor_qc/mark_duplicates_metrics
    tumor_insert_size_metrics:
        type: File[]
        outputSource: tumor_qc/insert_size_metrics
    tumor_alignment_summary_metrics:
        type: File[]
        outputSource: tumor_qc/alignment_summary_metrics
    tumor_hs_metrics:
        type: File[]
        outputSource: tumor_qc/hs_metrics
    tumor_per_target_coverage_metrics:
        type: File[]
        outputSource: tumor_qc/per_target_coverage_metrics
    tumor_per_target_hs_metrics:
        type: File[]
        outputSource: tumor_qc/per_target_hs_metrics
    tumor_per_base_coverage_metrics:
        type: File[]
        outputSource: tumor_qc/per_base_coverage_metrics
    tumor_per_base_hs_metrics:
        type: File[]
        outputSource: tumor_qc/per_base_hs_metrics
    tumor_summary_hs_metrics:
        type: File[]
        outputSource: tumor_qc/summary_hs_metrics
    tumor_flagstats:
        type: File[]
        outputSource: tumor_qc/flagstats
    tumor_verify_bam_id_metrics:
        type: File[]
        outputSource: tumor_qc/verify_bam_id_metrics
    tumor_verify_bam_id_depth:
        type: File[]
        outputSource: tumor_qc/verify_bam_id_depth
    
    mutect_unfiltered_vcf:
        type: File[]
        outputSource:  mutect/unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File[]
        outputSource: mutect/filtered_vcf
        secondaryFiles: [.tbi]
    mutect_nsamples_gnomad_pon_unfiltered:
        type: File
        outputSource: final_concatenated_vcf/final_concatenated_vcf
    mutect_nsamples_gnomad_pon_filtered:
        type: File
        outputSource: final_concatenated_vcf/final_concatenated_vcf
    mutect_merged_nsamples:
        type: File
        outputSource: combine_mutect_and_filter/merged_nsamples
        secondaryFiles: [.tbi]
    mutect_merged_nsamples_underN:
        type: File
        outputSource: combine_mutect_and_filter/merged_nsamples_underN
    mutect_merged_nsamples_overN:
        type: File
        outputSource: combine_mutect_and_filter/merged_nsamples_overN
        secondaryFiles: [.tbi]
    mutect_pon_total_counts:
        type: File
        outputSource: combine_mutect_and_filter/pon_total_counts
        secondaryFiles: [.tbi]
    mutect_final_concatenated:
        type: File
        outputSource: combine_mutect_and_filter/final_concatenated_vcf    
    
    varscan_unfiltered_vcf:
        type: File
        outputSource: detect_variants/varscan_unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: detect_variants/varscan_filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: detect_variants/pindel_unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: detect_variants/pindel_filtered_vcf
        secondaryFiles: [.tbi]
    docm_filtered_vcf:
        type: File
        outputSource: detect_variants/docm_filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: detect_variants/final_filtered_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: detect_variants/final_tsv
    vep_summary:
        type: File
        outputSource: detect_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_indel_bam_readcount_tsv
    

steps:
    tumor_cram_to_bam:
        scatter: [tumor_crams]
        run: ../subworkflows/cram_to_bam_and_index.cwl
        in:
            reference: reference
            tumor_crams: tumor_crams
        out: 
            [bam]
    bqsr:
        scatter: [bam]
        run: ../tools/bqsr.cwl
        in:
            reference: reference
            bam: tumor_cram_to_bam/bam
            intervals: bqsr_intervals
            known_sites: bqsr_known_sites
        out:
            [bqsr_table]
    apply_bqsr:
        scatter: [bam, bqsr_table]
        scatterMethod: dotproduct
        run: ../tools/apply_bqsr.cwl
        in:
            reference: reference
            bam: tumor_cram_to_bam/bam
            bqsr_table: bqsr/bqsr_table
            output_name: 
                valueFrom: "$(inputs.bam.nameroot).final"
        out:
            [bqsr_bam]  
    index_bam:
        scatter: [bam]
        run: ../tools/index_bam.cwl
        in:
            bam: apply_bqsr/bqsr_bam
        out:
            [indexed_bam] 
    tumor_qc
        scatter: [bam]
        run: ../subworkflows/qc_exome.cwl
        in:
            bam: index_bam/indexed_bam
            reference: reference
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
        out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    pad_target_intervals:
        run: ../tools/interval_list_expand.cwl
        in: 
            interval_list: target_intervals
            roi_padding: target_interval_padding
        out:
            [expanded_interval_list]
    ## scatter each caller
    mutect:
        scatter: [tumor_bam]
        run: ../subworkflows/mutect_normalize.cwl
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            scatter_count: scatter_count
            tumor_sample_name:
                valueFrom: "$(inputs.tumor_bam.nameroot)"
            filter_gnomAD_population_AF: filter_gnomAD_population_AF
        out:
            [unfiltered_vcf, filtered_vcf]
    mutect_on_target:
        scatter: [vcf]
        in:
            vcf = mutect/unfiltered_vcf
            targets = bed_targets
        out:
            [on_target_vcf]
    combine_mutect_and_filter:
        run: ../subworkflows/nsamples_and_PoN_and_filter.cwl
        in:
            reference: reference
            caller_vcfs: mutect_on_target/on_target_vcf
            nsamples: nsamples
            caller_prefix: "mutect2"
            normal_bams: pon_normal_bams
            pon_final_name: "mutect2.pon.total.counts"
        out:
            [final_concatenated_vcf, final_concatenated_filtered_vcf, merged_nsamples, merged_nsamples_underN, merged_nsamples_overN, pon_total_counts]
    mutect_annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: combine_mutect_and_filter/final_concatenated_filtered_vcf
            cache_dir: vep_cache_dir
            ensembl_assembly: vep_ensembl_assembly
            ensembl_version: vep_ensembl_version
            ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            coding_only: annotate_coding_only
            reference: reference
            pick: vep_pick
            custom_annotations: vep_custom_annotations
            plugins: vep_plugins
        out:
            [annotated_vcf, vep_summary]


    pindel:
        run: ../subworkflows/pindel.cwl
        scatter: [tumor_bam]
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            scatter_count: scatter_count
            insert_size: pindel_insert_size
            tumor_sample_name: 
                valueFrom: "$(inputs.tumor_bam.nameroot)"
            normal_sample_name: 
                valueFrom: "$(inputs.normal_bam.nameroot)"
        out:
            [unfiltered_vcf, filtered_vcf]
    pindel_on_target:
        scatter: [vcf]
        in:
            vcf = pindel/unfiltered_vcf
            targets = bed_targets
        out:
            [on_target_vcf]
    combine_pindel_and_filter:
