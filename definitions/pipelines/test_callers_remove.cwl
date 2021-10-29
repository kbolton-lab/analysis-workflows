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
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    tumor_name:
        type: string?
        default: 'tumor'
    tumor_sample_name:
        type: string
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    pon_normal_bams:
        type:
            type: array
            items: string
    scatter_count:
        type: int
        default: 20
    gatk_gnomad_af_only:
        type: File
        default:
            class: File
            path: "/storage1/fs1/bolton/Active/data/hg38/vcf/af-only-gnomad.biallelic.hg38.vcf.gz"
        secondaryFiles: [.tbi]
    filter_flag:
        type:
          - "null"
          - type: enum
            symbols: ["include", "exclude"]
        default: "include"
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
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
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    vep_plugins:
        type: string[]
        default: [Frameshift, Wildtype]
    af_threshold:
        type: float?
        default: 0.005
        doc: "for CH this is how low we want to call, used for vardict"
    bcbio_filter_string:
        type: string
        default: "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 45)))"
        doc: "http://bcb.io/2016/04/04/vardict-filtering/"
    pindel_insert_size:
        type: int
        default: 400
    ref_name:
        type: string?
        default: "GRCh38DH"
    ref_date:
        type: string?
        default: "20161216"
    pindel_min_supporting_reads:
        type: int?
        default: 3
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
    mutect_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    lofreq_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    vardict_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    pindel_pon2_file:
        type: File
        secondaryFiles: [.tbi]
    pon_pvalue:
        type: string?
        default: "4.098606e-08"
outputs:
    mutect_full:
        type: File
        outputSource: mutect/unfiltered_vcf
        doc: "full soft-filtered from fp_filter"
    mutect_pon_annotated_unfiltered_vcf:
        type: File
        outputSource: mutect_gnomad_pon_filters/processed_gnomAD_filtered_vcf
        doc: "gnomad filter only with fisher PoN p-value"
    mutect_pon_annotated_filtered_vcf:
        type: File
        outputSource: mutect_annotate_variants/annotated_vcf
        doc: "final annotated VCF with PoN fisher test hard filter"
    mutect_pon_total_counts:
        type: File
        outputSource: mutect_gnomad_pon_filters/pon_total_counts
        secondaryFiles: [.tbi]
        doc: "results from MSK pileup"

    # vardict_full:
    #     type: File
    #     outputSource: vardict/unfiltered_vcf
    #     doc: "full soft-filtered from fp_filter"
    # vardict_bcbio_filtered:
    #     type: File
    #     outputSource: vardict/bcbio_filtered_vcf
    #     secondaryFiles: [.tbi]
    #     doc: "bcbio soft filtered"
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
    #
    # pindel_full:
    #     type: File
    #     outputSource: pindel/unfiltered_vcf
    #     doc: "full soft-filtered from fp_filter"
    # pindel_pon_annotated_unfiltered_vcf:
    #     type: File
    #     outputSource: pindel_gnomad_pon_filters/processed_gnomAD_filtered_vcf
    #     doc: "gnomad filter only with fisher PoN p-value"
    # pindel_pon_annotated_filtered_vcf:
    #     type: File
    #     outputSource: pindel_pon2/annotated_vcf
    #     doc: "final annotated VCF with PoN fisher test hard filter"
    # pindel_pon_total_counts:
    #     type: File
    #     outputSource: pindel_gnomad_pon_filters/pon_total_counts
    #     secondaryFiles: [.tbi]
    #     doc: "results from MSK pileup"

    # lofreq_full:
    #     type: File
    #     outputSource: lofreq/unfiltered_vcf
    #     doc: "full soft-filtered from fp_filter"
    # lofreq_pon_annotated_unfiltered_vcf:
    #     type: File
    #     outputSource: lofreq_gnomad_pon_filters/processed_gnomAD_filtered_vcf
    #     doc: "gnomad filter only with fisher PoN p-value"
    # lofreq_pon_annotated_filtered_vcf:
    #     type: File
    #     outputSource: lofreq_pon2/annotated_vcf
    #     doc: "final annotated VCF with PoN fisher test hard filter"
    # lofreq_pon_total_counts:
    #     type: File
    #     outputSource: lofreq_gnomad_pon_filters/pon_total_counts
    #     secondaryFiles: [.tbi]
    #     doc: "results from MSK pileup"
    #
    # gnomad_exclude:
    #     type: File
    #     outputSource: get_gnomad_exclude/normalized_gnomad_exclude
    #     secondaryFiles: [.tbi]
steps:
    # bqsr:
    #     run: ../tools/bqsr.cwl
    #     in:
    #         reference: reference
    #         bam: tumor_bam
    #         intervals: bqsr_intervals
    #         known_sites: bqsr_known_sites
    #     out:
    #         [bqsr_table]
    # apply_bqsr:
    #     run: ../tools/apply_bqsr.cwl
    #     in:
    #         reference: reference
    #         bam: tumor_bam
    #         bqsr_table: bqsr/bqsr_table
    #         output_name:
    #             valueFrom: "$(inputs.bam.nameroot).final"
    #     out:
    #         [bqsr_bam]
    # index_bam:
    #     run: ../tools/index_bam.cwl
    #     in:
    #         bam: apply_bqsr/bqsr_bam
    #     out:
    #         [indexed_bam]
    # tumor_qc:
    #     run: ../subworkflows/qc_exome.cwl
    #     in:
    #         bam: index_bam/indexed_bam
    #         reference: reference
    #         bait_intervals: bait_intervals
    #         target_intervals: target_intervals
    #         per_base_intervals: per_base_intervals
    #         per_target_intervals: per_target_intervals
    #         summary_intervals: summary_intervals
    #         omni_vcf: omni_vcf
    #         picard_metric_accumulation_level: picard_metric_accumulation_level
    #         minimum_mapping_quality: qc_minimum_mapping_quality
    #         minimum_base_quality: qc_minimum_base_quality
    #     out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    # pad_target_intervals:
    #     run: ../tools/interval_list_expand.cwl
    #     in:
    #         interval_list: target_intervals
    #         roi_padding: target_interval_padding
    #     out:
    #         [expanded_interval_list]

    get_gnomad_exclude:
        run: ../subworkflows/get_gnomAD_filter.cwl
        in:
            reference: reference
            gnomad_AF_only: gatk_gnomad_af_only
            filter_flag: filter_flag
            output_type:
                default: "z"
            output_vcf_name:
                valueFrom: "gnomad.AF.exclude.vcf.gz"
        out: [gnomad_exclude, normalized_gnomad_exclude]
        doc: "this filter's the gnomAD_af_only file based on gnomAD POPAF threshold, it is what should be excluded if our calls have it since above threshold"

    # lofreq:
    #     run: ../subworkflows/lofreq.cwl
    #     in:
    #         reference: reference
    #         tumor_bam: tumor_bam
    #         interval_list: target_intervals
    #         scatter_count: scatter_count
    #         tumor_sample_name: tumor_sample_name
    #         min_var_freq: af_threshold
    #     out:
    #         [unfiltered_vcf, filtered_vcf]
    # lofreq_gnomad_pon_filters:
    #     run: ../subworkflows/gnomad_and_PoN_filter.cwl
    #     in:
    #         reference: reference
    #         caller_vcf: lofreq/unfiltered_vcf
    #         gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
    #         caller_prefix:
    #             source: tumor_sample_name
    #             valueFrom: "lofreq.$(self)"
    #         normal_bams: pon_normal_bams
    #         pon_final_name:
    #             source: tumor_sample_name
    #             valueFrom: "lofreq.$(self).pon.total.counts"
    #         pon_pvalue: pon_pvalue
    #     out:
    #         [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
    # lofreq_annotate_variants:
    #     run: ../tools/vep.cwl
    #     in:
    #         vcf: lofreq_gnomad_pon_filters/processed_filtered_vcf
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
    # lofreq_pon2:
    #     run: ../tools/pon2percent.cwl
    #     in:
    #         vcf: lofreq_annotate_variants/annotated_vcf
    #         vcf2PON: lofreq_pon2_file
    #         caller:
    #             valueFrom: "lofreq"
    #         sample_name: tumor_sample_name
    #     out:
    #         [annotated_vcf]

    # vardict:
    #     run: ../subworkflows/vardict_tumor_only.cwl
    #     in:
    #         reference: reference
    #         tumor_bam: tumor_bam
    #         interval_list: target_intervals # splits to bed
    #         scatter_count: scatter_count
    #         tumor_sample_name: tumor_sample_name
    #         af_threshold: af_threshold
    #         bcbio_filter_string: bcbio_filter_string
    #     out:
    #         [unfiltered_vcf, filtered_vcf, bcbio_filtered_vcf]
    # vardict_gnomad_pon_filters:
    #     run: ../subworkflows/gnomad_and_PoN_filter.cwl
    #     in:
    #         reference: reference
    #         caller_vcf: vardict/bcbio_filtered_vcf
    #         gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
    #         caller_prefix:
    #             source: tumor_sample_name
    #             valueFrom: "vardict.$(self)"
    #         normal_bams: pon_normal_bams
    #         pon_final_name:
    #             source: tumor_sample_name
    #             valueFrom: "vardict.$(self).pon.total.counts"
    #         pon_pvalue: pon_pvalue
    #     out:
    #         [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
    # vardict_annotate_variants:
    #     run: ../tools/vep.cwl
    #     in:
    #         vcf: vardict_gnomad_pon_filters/processed_filtered_vcf
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
    # vardict_pon2:
    #     run: ../tools/pon2percent.cwl
    #     in:
    #         vcf: vardict_annotate_variants/annotated_vcf
    #         vcf2PON: vardict_pon2_file
    #         caller:
    #             valueFrom: "vardict"
    #         sample_name: tumor_sample_name
    #     out:
    #         [annotated_vcf]
    #
    # pindel:
    #     run: ../subworkflows/pindel_tumor_only.cwl
    #     in:
    #         reference: reference
    #         tumor_bam: tumor_bam
    #         interval_list: target_intervals
    #         scatter_count: scatter_count
    #         insert_size: pindel_insert_size
    #         tumor_sample_name: tumor_sample_name
    #         ref_name: ref_name
    #         ref_date: ref_date
    #         pindel_min_supporting_reads: pindel_min_supporting_reads
    #         min_var_freq: af_threshold
    #     out:
    #         [unfiltered_vcf, filtered_vcf]
    # pindel_gnomad_pon_filters:
    #     run: ../subworkflows/gnomad_and_PoN_filter.cwl
    #     in:
    #         reference: reference
    #         caller_vcf: pindel/unfiltered_vcf
    #         gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
    #         caller_prefix:
    #             source: tumor_sample_name
    #             valueFrom: "pindel.$(self)"
    #         normal_bams: pon_normal_bams
    #         pon_final_name:
    #             source: tumor_sample_name
    #             valueFrom: "pindel.$(self).pon.total.counts"
    #         pon_pvalue: pon_pvalue
    #     out:
    #         [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
    # pindel_annotate_variants:
    #     run: ../tools/vep.cwl
    #     in:
    #         vcf: pindel_gnomad_pon_filters/processed_filtered_vcf
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
    # pindel_pon2:
    #     run: ../tools/pon2percent.cwl
    #     in:
    #         vcf: pindel_annotate_variants/annotated_vcf
    #         vcf2PON: pindel_pon2_file
    #         caller:
    #             valueFrom: "pindel"
    #         sample_name: tumor_sample_name
    #     out:
    #         [annotated_vcf]
    #
    mutect:
        run: ../subworkflows/mutect_no_scatter.cwl
        in:
            reference: reference
            tumor_bam: index_bam/indexed_bam
            # interval_list: pad_target_intervals/expanded_interval_list
            interval_list: target_intervals
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    mutect_gnomad_pon_filters:
        run: ../subworkflows/gnomad_and_PoN_filter.cwl
        in:
            reference: reference
            caller_vcf: mutect/unfiltered_vcf
            gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
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
    mutect_annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: mutect_gnomad_pon_filters/processed_filtered_vcf
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
