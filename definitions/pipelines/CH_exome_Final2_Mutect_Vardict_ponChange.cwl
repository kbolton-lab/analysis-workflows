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
    pon_normal_bams:
        type: 
            type: array
            items: string
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
          In general for a WES exome reagent bait_intervals and target_intervals are the same. UKBB is xgen_plus_spikein.GRCh38.intervals
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
    scatter_count:
        type: int
        default: 20
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
    impact_annotation:
        type: File
    topmed_annotation:
        type: File
    cosmic_annotation:
        type: File
    tsg_annotation:
        type: File
    oncoKB_annotation:
        type: File
    pd_table_annotation:
        type: File
    panmyeloid_annotation:
        type: File
    blacklist_annotation:
        type: File
    segemental_duplications_annotation:
        type: File
    simple_repeats_annotation:
        type: File
    repeat_masker_annotation:
        type: File
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
        default: 0.005
    gatk_gnomad_af_only:
        type: File?
        default: 
            class: File
            path: "/cache2/fs1/bolton/Active/data/hg38/vcf/af-only-gnomad.biallelic.hg38.vcf.gz"
        secondaryFiles: [.tbi]
    gatk_gnomad_af_only_hc_0_005:
        type: File?
        default: 
            class: File
            path: "/cache2/fs1/bolton/Active/data/hg38/vcf/af-only-gnomad.biallelic.above.005.hg38.vcf.gz"
        secondaryFiles: [.tbi]
    filter_flag:
        type:
          - "null" 
          - type: enum
            symbols: ["include", "exclude"]
        default: "include"
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    cle_vcf_filter:
        type: boolean
        default: false
    variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
        doc: "The names of one or more standard VCF fields or INFO fields to include in the output table. NOT USING"
    variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
        doc: "The name of a genotype field to include in the output table. NOT USING"
    vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
        doc: "VEP fields in final output. NOT USING"
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "An optional VCF with variants that will be flagged as 'VALIDATED' if found in this pipeline's main output VCF"
    vep_plugins:
        type: string[]
        default: [Frameshift, Wildtype]
    af_threshold: 
        type: float?
        default: 0.005
        doc: "for CH this is how low we want to call, used for vardict"
    bcbio_filter_string:
        type: string
        default: "((FMT/AF * FMT/DP < 6) && ((FMT/MQ < 55.0 && FMT/NM > 1.0) || (FMT/MQ < 60.0 && FMT/NM > 2.0) || (FMT/DP < 10) || (FMT/QUAL < 45)))"
        doc: "http://bcb.io/2016/04/04/vardict-filtering/"
    mutect_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    varscan_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    vardict_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    pindel_pon2_file:
        type: File
        secondaryFiles: [.tbi]


outputs:
    tumor_alignment_summary_metrics:
        type: File
        outputSource: tumor_qc/alignment_summary_metrics
    tumor_hs_metrics:
        type: File
        outputSource: tumor_qc/hs_metrics
    
    mutect_full:
        type: File
        outputSource: mutect/unfiltered_vcf
        doc: "FilterMutectCalls soft filtered"
    mutect_annotated_vcf:
        type: File
        outputSource: mutect_pon2/annotated_vcf
        doc: "VEP, mutect soft filtered, gnomad filter"

    # mutect_pon_annotated_unfiltered_vcf:
    #     type: File
    #     outputSource: mutect_gnomad_pon_filters/processed_gnomAD_filtered_vcf
    #     doc: "gnomad filter only with fisher PoN p-value"
    # mutect_pon_annotated_filtered_vcf:
    #     type: File
    #     outputSource: mutect_pon2/annotated_vcf
    #     doc: "final annotated VCF with PoN fisher test hard filter"
    # mutect_pon_total_counts:
    #     type: File
    #     outputSource: mutect_gnomad_pon_filters/pon_total_counts
    #     secondaryFiles: [.tbi]
    #     doc: "results from MSK pileup"


    vardict_bcbio_filtered:
        type: File
        outputSource: vardict/bcbio_filtered_vcf
        secondaryFiles: [.tbi]
        doc: "bcbio soft filtered"
    vardict_annotated_vcf:
        type: File
        outputSource: vardict_pon2/annotated_vcf
        doc: "VEP, bcbio soft filtered, gnomad filter"
    
    # vardict_pon_annotated_unfiltered_vcf:
    #     type: File
    #     outputSource: vardict_gnomad_pon_filters/processed_gnomAD_filtered_vcf
    #     doc: "gnomad filter only with fisher PoN p-value"
    # vardict_pon_annotated_filtered_vcf:
    #     type: File
    #     outputSource: vardict_pon2/annotated_vcf
    #     doc: "final annotated VCF with PoN fisher test hard filter"
    # vardict_pon_total_counts:
    #     type: File
    #     outputSource: vardict_gnomad_pon_filters/pon_total_counts
    #     secondaryFiles: [.tbi]
    #     doc: "results from MSK pileup"

    final_tsv: 
        type: File
        outputSource: final_annotation/final_tsv
    column_check:
        type: File
        outputSource: final_annotation/column_check

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
    bqsr:
        run: ../tools/bqsr.cwl
        in:
            reference: reference
            bam: tumor_cram_to_bam/bam
            intervals: bqsr_intervals
            known_sites: bqsr_known_sites
        out:
            [bqsr_table]
    apply_bqsr:
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
        run: ../tools/index_bam.cwl
        in:
            bam: apply_bqsr/bqsr_bam
        out:
            [indexed_bam] 
    tumor_qc:
        run: ../subworkflows/qc_exome_UKBB.cwl
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
        # out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
        out: [alignment_summary_metrics, hs_metrics]

    mutect:
        run: ../subworkflows/mutect_normalize_no_fp.cwl
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            normal_bam: normal_bam
            interval_list: target_intervals
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
        out:
            [unfiltered_vcf]
    mutect_extract_tumor:
        run: ../tools/bcftools_extract_tumor.cwl
        in:
            vcf: mutect/unfiltered_vcf
            output_type:
                default: "z"
            tumor_sample_name: tumor_sample_name
            output_vcf_name: 
                valueFrom: $(inputs.tumor_sample_name).vcf.gz
        out:
            [tumor_only_vcf]
    mutect_extract_tumor_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: mutect_extract_tumor/tumor_only_vcf
        out:
            [indexed_vcf]   
    mutect_gnomad_pon_filters:
        run: ../subworkflows/gnomad_and_PoN_filter.cwl
        in:
            reference: reference
            caller_vcf: mutect_extract_tumor_index/indexed_vcf
            gnomAD_exclude_vcf: gatk_gnomad_af_only_hc_0_005
            caller_prefix: 
                source: tumor_sample_name
                valueFrom: "mutect.$(self)"
            normal_bams: pon_normal_bams
            pon_final_name: 
                source: tumor_sample_name
                valueFrom: "mutect.$(self).pon.total.counts"
        out:
            [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
        doc: "processed_filtered_vcf is gnomAD and PoN filtered"
    # Mutect_isec_complement_gnomAD:
    #     run: ../tools/bcftools_isec_complement.cwl
    #     in:
    #         vcf: mutect_extract_tumor_index/indexed_vcf
    #         exclude_vcf: gatk_gnomad_af_only_hc_0_005
    #         output_vcf_name: 
    #             valueFrom: "mutect.gnomAD_AF_filter.vcf.gz"
    #         output_type: 
    #             default: "z"
    #     doc: "both input files should be bgzipped and indexed"
    #     out:
    #         [complement_vcf]
    # mutect_add_sample_info:
    #     run: ../sample_info_tag.cwl
    #     in: 
    #         Mutect_isec_complement_gnomAD/complement_vcf
    #         output_vcf_name: 
    #             valueFrom: "mutect.gnomAD_AF_filter.sample.vcf.gz"
    #         output_type: 
    #             default: "z"
    #          out:
    #             [sample_vcf]
    mutect_annotate_variants:
        run: ../tools/vep.cwl
        in:
            #vcf: mutect_add_sample_info/sample_vcf
            vcf: mutect_gnomad_pon_filters/processed_gnomAD_filtered_vcf
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
    mutect_pon2:
        run: ../tools/pon2percent.cwl
        in:
            vcf: mutect_annotate_variants/annotated_vcf
            vcf2PON: mutect_pon2_file
            caller:  
                valueFrom: "mutect"
            sample_name: tumor_sample_name
        out:
            [annotated_vcf]

    vardict:
        run: ../subworkflows/vardict_no_fp.cwl
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            normal_bam: normal_bam
            interval_list: target_intervals # splits to bed
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
            af_threshold: af_threshold
            normal_sample_name: normal_sample_name
            bcbio_filter_string: bcbio_filter_string
        out:
            [bcbio_filtered_vcf]
    Vardict_isec_complement_gnomAD:
        run: ../tools/bcftools_isec_complement.cwl
        in:
            vcf: vardict/bcbio_filtered_vcf
            exclude_vcf: gatk_gnomad_af_only_hc_0_005
            output_vcf_name: 
                valueFrom: "vardict.gnomAD_AF_filter.vcf.gz"
            output_type: 
                default: "z"
        doc: "both input files should be bgzipped and indexed"
        out:
            [complement_vcf]
    Vardict_gnomAD_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: Vardict_isec_complement_gnomAD/complement_vcf
        out:
            [indexed_vcf]   
    vardict_intersect_mutect:
        run: ../tools/bcftools_isec.cwl
        in:
            vcf1: Vardict_gnomAD_index/indexed_vcf
            vcf2: mutect_pon2/annotated_vcf
            output_type: 
                default: "z"
            output_vcf_name: 
                valueFrom: "vardict.intersect.MUTECT2.vcf.gz"
        out:
            [overlap_vcf]
    # vardict_annotate_variants:
    #     run: ../tools/vep.cwl
    #     in:
    #         vcf: Vardict_isec_complement_gnomAD/complement_vcf
    #         cache_dir: vep_cache_dir
    #         ensembl_assembly: vep_ensembl_assembly
    #         ensembl_version: vep_ensembl_version
    #         ensembl_species: vep_ensembl_species
    #         synonyms_file: synonyms_file
    #         coding_only: annotate_coding_only
    #         reference: reference
    #         pick: vep_pick
    #         custom_annotations: vep_custom_annotations
    #         plugins: vep_plugins
    #     out:
    #         [annotated_vcf, vep_summary]
    vardict_pon2:
        run: ../tools/pon2percent.cwl
        in:
            vcf: vardict_intersect_mutect/overlap_vcf
            vcf2PON: vardict_pon2_file
            caller:  
                valueFrom: "vardict"
            sample_name: tumor_sample_name
        out:
            [annotated_vcf]

    final_annotation:
        run: ../tools/annotate_CH_pd_Mutect_Vardict_ponChange.cwl
        in:
            mutect_vcf: mutect_pon2/annotated_vcf
            vardict_vcf: vardict_pon2/annotated_vcf
            sample_name: tumor_sample_name
            impact_annotation: impact_annotation
            topmed_annotation: topmed_annotation
            cosmic_annotation: cosmic_annotation
            tsg_annotation: tsg_annotation
            oncoKB_annotation: oncoKB_annotation
            pd_table_annotation: pd_table_annotation
            panmyeloid_annotation: panmyeloid_annotation
            blacklist_annotation: blacklist_annotation
            segemental_duplications_annotation: segemental_duplications_annotation
            simple_repeats_annotation: simple_repeats_annotation
            repeat_masker_annotation: repeat_masker_annotation
            intervals: target_intervals
        out:
            [final_tsv, column_check]
    