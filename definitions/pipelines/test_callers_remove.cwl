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
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
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
            path: "/storage1/fs1/bolton/Active/data/hg19/vcf/af-only-gnomad.raw.sites.vcf.gz"
        secondaryFiles: [.tbi]
    filter_flag:
        type:
          - "null"
          - type: enum
            symbols: ["include", "exclude"]
        default: "include"
        doc: "we default the gnomad to include > than threshold because this is smaller file for isec complement"
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
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

    vardict_full:
        type: File
        outputSource: vardict/unfiltered_vcf
        doc: "full soft-filtered from fp_filter"
    vardict_bcbio_filtered:
        type: File
        outputSource: vardict/bcbio_filtered_vcf
        secondaryFiles: [.tbi]
        doc: "bcbio soft filtered"
    vardict_pon_annotated_unfiltered_vcf:
        type: File
        outputSource: vardict_gnomad_pon_filters/processed_gnomAD_filtered_vcf
        doc: "gnomad filter only with fisher PoN p-value"
    vardict_pon_annotated_filtered_vcf:
        type: File
        outputSource: vardict_annotate_variants/annotated_vcf
        doc: "final annotated VCF with PoN fisher test hard filter"
    vardict_pon_total_counts:
        type: File
        outputSource: vardict_gnomad_pon_filters/pon_total_counts
        secondaryFiles: [.tbi]
        doc: "results from MSK pileup"

    pindel_full:
        type: File
        outputSource: pindel/unfiltered_vcf
        doc: "full soft-filtered from fp_filter"
    pindel_pon_annotated_unfiltered_vcf:
        type: File
        outputSource: pindel_gnomad_pon_filters/processed_gnomAD_filtered_vcf
        doc: "gnomad filter only with fisher PoN p-value"
    pindel_pon_annotated_filtered_vcf:
        type: File
        outputSource: pindel_annotate_variants/annotated_vcf
        doc: "final annotated VCF with PoN fisher test hard filter"
    pindel_pon_total_counts:
        type: File
        outputSource: pindel_gnomad_pon_filters/pon_total_counts
        secondaryFiles: [.tbi]
        doc: "results from MSK pileup"

    gnomad_exclude:
        type: File
        outputSource: get_gnomad_exclude/normalized_gnomad_exclude
        secondaryFiles: [.tbi]
steps:
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

    pindel:
        run: ../subworkflows/pindel_tumor_only.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            interval_list: target_intervals
            scatter_count: scatter_count
            insert_size: pindel_insert_size
            tumor_sample_name: tumor_sample_name
            ref_name: ref_name
            ref_date: ref_date
            pindel_min_supporting_reads: pindel_min_supporting_reads
        out:
            [unfiltered_vcf, filtered_vcf]
    pindel_gnomad_pon_filters:
        run: ../subworkflows/gnomad_and_PoN_filter.cwl
        in:
            reference: reference
            caller_vcf: pindel/unfiltered_vcf
            gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
            caller_prefix:
                source: tumor_sample_name
                valueFrom: "pindel.$(self)"
            normal_bams: pon_normal_bams
            pon_final_name:
                source: tumor_sample_name
                valueFrom: "pindel.$(self).pon.total.counts"
        out:
            [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
    pindel_annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: pindel_gnomad_pon_filters/processed_filtered_vcf
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
    vardict:
        run: ../subworkflows/vardict_tumor_only.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            interval_list: target_intervals # splits to bed
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
            af_threshold: af_threshold
            bcbio_filter_string: bcbio_filter_string
        out:
            [unfiltered_vcf, filtered_vcf, bcbio_filtered_vcf]
    vardict_gnomad_pon_filters:
        run: ../subworkflows/gnomad_and_PoN_filter.cwl
        in:
            reference: reference
            caller_vcf: vardict/bcbio_filtered_vcf
            gnomAD_exclude_vcf: get_gnomad_exclude/normalized_gnomad_exclude
            caller_prefix:
                source: tumor_sample_name
                valueFrom: "vardict.$(self)"
            normal_bams: pon_normal_bams
            pon_final_name:
                source: tumor_sample_name
                valueFrom: "vardict.$(self).pon.total.counts"
        out:
            [processed_gnomAD_filtered_vcf, processed_filtered_vcf, pon_total_counts]
    vardict_annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: vardict_gnomad_pon_filters/processed_filtered_vcf
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

    mutect:
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
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
