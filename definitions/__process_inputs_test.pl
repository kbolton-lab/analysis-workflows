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
            mutect_pon2 => { input_type => 'File' },
            vardict_pon2 => { input_type => 'File' },
            varscan_pon2 => { input_type => 'File' },
            pindel_pon2 => { input_type => 'File' },
            target_intervals => { input_type => 'File' },
            bait_intervals => { input_type => 'File' },
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
$inputs->{scatter_count} = 50;
$inputs->{bqsr_intervals} = [map 'chr'.$_, (1..22,X,Y,M)];

# pon
my ($pon_file, @extra) = grep { $_->name eq 'pon_normal_bams' } @inputs;
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


## What is this section? #################################################################
# my ($tumor_input, @extra) = grep { $_->name eq 'tumor_sample' } @inputs;
# if (@extra) {
#     die 'multiple inputs found for tumor sample';
# }
# my $tumor_sample = $tumor_input->value_class_name->get($tumor_input->value_id)
#     or die 'no tumor found for input';

# my ($normal_input, @extra) = grep { $_->name eq 'normal_sample' } @inputs;
# if (@extra) {
#     die 'multiple inputs found for normal sample';
# }
# my $normal_sample = $normal_input->value_class_name->get($normal_input->value_id)
#     or die 'no normal found for input';
##########################################################################################

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
my @ids = map $_->value_id, @instrument_data_inputs;

#my $target_region_set_name;
my $multiple_trsn = 0;
$inputs->{tumor_sample_name} = $build->subject->name;

# cat /gscmnt/gc2560/core/processing-profile/cwl-pipeline/f6faba8d6c234cac9bcd58e60ebe7579/process_inputs.pl
# cat /gscmnt/gc2560/core/processing-profile/cwl-pipeline/f5eb943e8a7b4fcbb7a4c26e4b98e1a3/process_inputs.pl
my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
# my @ids = map $_->value_id, @instrument_data_inputs; # this is never used even in somatic_exome process_inputs.pl

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
## only 1 unmatched normal
my ($normal_bam, @extra) = grep { $_->name eq 'normal_bam' } @inputs;
if (@extra) {
    die 'multiple inputs found for normal bam';
$inputs->{normal_bam} = {class => 'File', path => $normal_bam->value_id}; 
my ($normal_sample_name, @extra) = grep { $_->name eq 'normal_sample_name' } @inputs;
if (@extra) {
    die 'multiple inputs found for normal bam';
$inputs->{normal_sample_name} = $normal_sample_name->value_id;



$inputs->{per_base_intervals} = [
    { label => 'clinvar', file => '/gscmnt/gc2560/core/model_data/interval-list/01f4fae3699646c3af2fa47853da7a8c/06a82ecf9c434b7ab03d82e59eaa28c8.interval_list' },
];
$inputs->{per_target_intervals} = [
    { label => 'acmg_genes', file => '/gscmnt/gc2560/core/model_data/interval-list/db8c25932fd94d2a8a073a2e20449878/a35b64d628b94df194040032d53b5616.interval_list' },
];
$inputs->{summary_intervals} = [];

$inputs->{variants_to_table_fields} = [qw(CHROM POS REF ALT set)];
$inputs->{variants_to_table_genotype_fields} = [qw(GT AD AF DP)];
$inputs->{vep_to_table_fields} = [qw(Consequence SYMBOL Feature_type Feature HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons HGNC_ID Existing_variation gnomADe_AF CLIN_SIG SOMATIC PHENO clinvar_CLINSIGN clinvar_PHENOTYPE clinvar_SCORE clinvar_RCVACC clinvar_TESTEDINGTR clinvar_PHENOTYPELIST clinvar_NUMSUBMIT clinvar_GUIDELINES)];

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