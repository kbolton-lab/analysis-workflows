#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Spec;
use YAML::XS;

my $build_id = $ARGV[0]
    or die 'no build id';

my $build = Genome::Model::Build->get($build_id)
    or die 'no build for id';

{
    package InputProcessor;
    class InputProcessor {
        is => 'Genome::Model::Build::CwlPipeline::InputProcessor',
        has_simple_input => [
            reference => { input_type => 'Text' },
            mutect_pon2_file => { input_type => 'File' },
            vardict_pon2_file => { input_type => 'File' },
            varscan_pon2_file => { input_type => 'File' },
            pindel_pon2_file => { input_type => 'File' },
            target_intervals => { input_type => 'File' },
            bait_intervals => { input_type => 'File' },
            gatk_gnomad_af_only_hc_0_005 => { input_type => 'File' },
            impact_annotation => { input_type => 'File' },
            topmed_annotation => { input_type => 'File' },
            cosmic_annotation => { input_type => 'File' },
            tsg_annotation => { input_type => 'File' },
            oncoKB_annotation => { input_type => 'File' },
            pd_table_annotation => { input_type => 'File' },
            panmyeloid_annotation => { input_type => 'File' },
            blacklist_annotation => { input_type => 'File' },
            segemental_duplications_annotation => { input_type => 'File' },
            simple_repeats_annotation => { input_type => 'File' },
            repeat_masker_annotation => { input_type => 'File' },
            picard_metric_accumulation_level => { input_type => 'Text' },
            vep_cache_dir => { input_type => 'Text' },
            mills => { input_type => 'File' },
            known_indels => { input_type => 'File' },
            dbsnp_vcf => { input_type => 'File' },
            docm_vcf => { input_type => 'File' },
            synonyms_file => { input_type => 'File' },
            omni_vcf => { input_type => 'File' },
            cosmic_vcf => { input_type => 'File' },
            panel_of_normals_vcf => { input_type => 'File' },
            vep_assembly => { input_type => 'Text' },
            manta_call_regions => { input_type => 'File' },
            manta_non_wgs => { input_type => 'Text' },
            manta_output_contigs => { input_type => 'Text' },
            somalier_vcf => { input_type => 'File' },
            vep_ensembl_assembly => { input_type => 'Text' },
            vep_ensembl_version => { input_type => 'Text' },
            vep_ensembl_species => { input_type => 'Text' }
        ],
    };
}

my @inputs = $build->inputs;
my $input_processor = InputProcessor->get($build->id);
my $inputs = $input_processor->simple_inputs;
$inputs->{scatter_count} = 20;
$inputs->{bqsr_intervals} = [map 'chr'.$_, (1..22,"X","Y","M")];

# pon
my ($pon_file, @extra) = grep { $_->name eq 'pon_normal_bams_bulk' } @inputs;
if (@extra) {
    die 'multiple inputs found for PoN file';
}
my @pon_lines = Genome::Sys->read_file($pon_file->value_id); chomp @pon_lines;
my @pon_bams;
foreach (@pon_lines) {
    #push @pon_bams, { class => 'File', path => $_ };
    ## make string instead to list of /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/young/<sample>.bam
    push @pon_bams, $_
}
$inputs->{pon_normal_bams} = \@pon_bams;

# bqsr_known_sites
# my ($bqsr_file, @extra) = grep { $_->name eq 'bqsr_known_sites' } @inputs;
# if (@extra) {
#     die 'multiple inputs found for bqsr_known_sites file';
# }
# my @bqsr_lines = Genome::Sys->read_file($bqsr_file->value_id); chomp @bqsr_lines;
# my @bqsr_known_sites;
# foreach (@bqsr_lines) {
#     push @bqsr_known_sites, { class => 'File', path => $_ };
# }
# $inputs->{bqsr_known_sites} = \@bqsr_known_sites;
my @bqsr_known_sites;
for my $input_name (qw(dbsnp_vcf known_indels mills)) {
    my ($file_input, @extra) = grep { $_->name eq $input_name } @inputs;
    if (@extra) { die 'multiple inputs found for ' . $input_name; }
    push @bqsr_known_sites, { class => 'File', path => $file_input->value_id };
}

$inputs->{bqsr_known_sites} = \@bqsr_known_sites;

# my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
# my @ids = map $_->value_id, @instrument_data_inputs;

#my $target_region_set_name;
my $multiple_trsn = 0;
$inputs->{tumor_sample_name} = $build->subject->name;

# cat /gscmnt/gc2560/core/processing-profile/cwl-pipeline/f6faba8d6c234cac9bcd58e60ebe7579/process_inputs.pl
# cat /gscmnt/gc2560/core/processing-profile/cwl-pipeline/f5eb943e8a7b4fcbb7a4c26e4b98e1a3/process_inputs.pl
# my @ids = map $_->value_id, @instrument_data_inputs; # this is never used even in somatic_exome process_inputs.pl
my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    if (my @crams = sort glob(File::Spec->join($id->disk_allocation->absolute_path, '*.cram'))) {
        # there is only one cram per sample ever
       $inputs->{tumor_cram} = { class => 'File', path => $crams[0] };
    } else {
        die 'No CRAM found for instrument data ' . $id->id;
    } 
}
######################################################
## only 1 unmatched normal input for ever sample CRAM
my ($normal_bam, @extra) = grep { $_->name eq 'normal_bam' } @inputs;
if (@extra) {
    die 'multiple inputs found for normal bam';
}
$inputs->{normal_bam} = {class => 'File', path => $normal_bam->value_id}; 
my ($normal_sample_name, @extra) = grep { $_->name eq 'normal_sample_name' } @inputs;
if (@extra) {
    die 'multiple inputs found for normal bam';
}
$inputs->{normal_sample_name} = $normal_sample_name->value_id;


$inputs->{per_base_intervals} = [
    { label => 'UKBB_WES', file => '/cache2/fs1/bolton/Active/projects/mocha/UKBB/exome/xgen_plus_spikein.GRCh38.interval_list' },
];
$inputs->{per_target_intervals} = [
    { label => 'UKBB_WES', file => '/cache2/fs1/bolton/Active/projects/mocha/UKBB/exome/xgen_plus_spikein.GRCh38.interval_list' },
];
$inputs->{summary_intervals} = [];

$inputs->{variants_to_table_fields} = [qw(CHROM POS REF ALT set)];
$inputs->{variants_to_table_genotype_fields} = [qw(GT AD AF DP)];
## adding gnomADg_AF for gnomAD 3.1 genomes
$inputs->{vep_to_table_fields} = [qw(Consequence SYMBOL Feature_type Feature HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons HGNC_ID Existing_variation gnomADe_AF gnomADg_AF CLIN_SIG SOMATIC PHENO clinvar_CLINSIGN clinvar_PHENOTYPE clinvar_SCORE clinvar_RCVACC clinvar_TESTEDINGTR clinvar_PHENOTYPELIST clinvar_NUMSUBMIT clinvar_GUIDELINES)];

my @vep_custom_annotations;

my ($gnomad_file, @extra) = grep { $_->name eq 'custom_gnomad_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_gnomad_vcf';
}
if($gnomad_file) {
    my $vcf_path = $gnomad_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'gnomADe'; # if changed should update the `vep_to_table_fields` input below
    $annotation_info->{gnomad_filter} = 'true';
    $annotation_info->{check_existing} = 'true';
    $annotation_info->{vcf_fields} = ['AF','AF_AFR','AF_AMR','AF_ASJ','AF_EAS','AF_FIN','AF_NFE','AF_OTH','AF_SAS'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}

my ($gnomad_file, @extra) = grep { $_->name eq 'custom_gnomadV3_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_gnomadV3_vcf';
}
if($gnomad3_file) {
    my $vcf_path = $gnomad3_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'gnomADg'; # if changed should update the `vep_to_table_fields` input below
    $annotation_info->{gnomad_filter} = 'true';
    $annotation_info->{check_existing} = 'true';
    $annotation_info->{vcf_fields} = ['AF','AF_ami','AF_oth','AF_afr','AF_sas','AF_asj','AF_fin','AF_amr','AF_nfe','AF_eas'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}

my ($clinvar_file, @extra) = grep { $_->name eq 'custom_clinvar_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_clinvar_vcf';
}
if($clinvar_file) {
    my $vcf_path = $clinvar_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'clinvar';
    $annotation_info->{gnomad_filter} = 'false';
    $annotation_info->{check_existing} = 'false';
    $annotation_info->{vcf_fields} = ['CLINSIGN','PHENOTYPE','SCORE','RCVACC','TESTEDINGTR','PHENOTYPELIST','NUMSUBMIT','GUIDELINES'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}
$inputs->{vep_custom_annotations} = \@vep_custom_annotations;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, $inputs);