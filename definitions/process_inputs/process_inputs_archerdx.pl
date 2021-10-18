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
            UMI_paired => { input_type => 'Number' },
            umi_length => { input_type => 'Number' },
            max_read_error_rate => { input_type => 'Number' },
            max_base_error_rate => { input_type => 'Number' },
            min_base_quality => { input_type => 'Number' },
            max_no_call_fraction => { input_type => 'Number' },
            af_threshold => { input_type => 'Number' },
            pon_pvalue => { input_type => 'Text' },
            tumor_name => { input_type => 'Text' },
            mutect_pon2_file => { input_type => 'File' },
            vardict_pon2_file => { input_type => 'File' },
            varscan_pon2_file => { input_type => 'File' },
            pindel_pon2_file => { input_type => 'File' },
            target_intervals => { input_type => 'File' },
            target_interval_padding => { input_type => 'Number' },
            bait_intervals => { input_type => 'File' },
            interval_list => { input_type => 'File' },
            qc_minimum_mapping_quality => { input_type => 'Number' },
            qc_minimum_base_quality => { input_type => 'Number' },
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
            mills => { input_type => 'File' },
            known_indels => { input_type => 'File' },
            dbsnp_vcf => { input_type => 'File' },
            docm_vcf => { input_type => 'File' },
            synonyms_file => { input_type => 'File' },
            omni_vcf => { input_type => 'File' },
            cosmic_vcf => { input_type => 'File' },
            panel_of_normals_vcf => { input_type => 'File' },
            somalier_vcf => { input_type => 'File' },
            vep_cache_dir => { input_type => 'Text' },
            vep_ensembl_assembly => { input_type => 'Text' },
            vep_ensembl_version => { input_type => 'Text' },
            vep_ensembl_species => { input_type => 'Text' },
            validated_variants => { input_type => 'File' },
            mutect_max_alt_allele_in_normal_fraction => { input_type =>'Number' },
            mutect_max_alt_alleles_in_normal_count => { input_type => 'Number' },
            pindel_insert_size => { inputy_type => 'Number' },
            ref_name => { input_type => 'Text' },
            ref_date => { input_type => 'Text' },
            pindel_min_supporting_reads => { input_type => 'Number' },
            filter_gnomADe_maximum_population_allele_frequency => { input_type => 'Number' }
        ],
    };
}

my @inputs = $build->inputs;
my $input_processor = InputProcessor->get($build->id);
my $inputs = $input_processor->simple_inputs;
$inputs->{scatter_count} = 50;
$inputs->{bqsr_intervals} = [map 'chr'.$_, (1..22,"X","Y","M")];

# Sequence
my @sequence;
my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
my @ids = map $_->value_id, @instrument_data_inputs;

my $target_region_set_name;
my $multiple_trsn = 0;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'Archer';
    my $rgid = $id->id;

    my $sequence = {};
    $sequence->{readgroup} = join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");

    if (my $bam_path = $id->bam_path) {
        $sequence->{sequence}{bam} = {class => 'File', path => $bam_path};
    } elsif (my @fastqs = sort glob(File::Spec->join($id->disk_allocation->absolute_path, '*.fastq.gz'))) {
        unless (@fastqs == 2) {
            die "expected two fastqs but got " . scalar(@fastqs);
        }

        $sequence->{sequence}{fastq1} = { class => 'File', path => $fastqs[0] };
        $sequence->{sequence}{fastq2} = { class => 'File', path => $fastqs[1] };
    } else {
        die 'No FASTQs or BAM found for instrument data ' . $id->id;
    }

    push @sequence, $sequence;

    unless ($target_region_set_name) {
        $target_region_set_name = $id->target_region_set_name;
    } else {
        if ($id->target_region_set_name ne $target_region_set_name) {
            $multiple_trsn = 1;
        }
    }
}

$inputs->{sequence} = \@sequence;

# Read Structure
my @read_structure_inputs = Genome::Model::Build::Input->get(build_id => $build->id, 'name LIKE' => 'read_structure_R%', -order_by => 'name');
my @read_structure = map $_->value_id, @read_structure_inputs;

$inputs->{read_structure} = \@read_structure;

# Min Reads
my @min_reads_input = grep { $_->name eq 'min_reads' } @inputs;
my @min_reads = map { 0 + $_->value_id } @min_reads_input;

$inputs->{min_reads} = \@min_reads;

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

$inputs->{per_base_intervals} = [
    { label => 'Archer_Panel', file => '/storage1/fs1/bolton/Active/data/hg38/bed/archer_panel.hg38.interval_list' },
];
$inputs->{per_target_intervals} = [
    { label => 'Archer_Panel', file => '/storage1/fs1/bolton/Active/data/hg38/bed/archer_panel.hg38.interval_list' },
];
$inputs->{summary_intervals} = [];

#$inputs->{variants_to_table_fields} = [qw(CHROM POS REF ALT set)];
#$inputs->{variants_to_table_genotype_fields} = [qw(GT AD AF DP)];
#$inputs->{vep_to_table_fields} = [qw(Consequence SYMBOL Feature_type Feature HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons HGNC_ID Existing_variation gnomADe_AF CLIN_SIG SOMATIC PHENO clinvar_CLINSIGN clinvar_PHENOTYPE clinvar_SCORE clinvar_RCVACC clinvar_TESTEDINGTR clinvar_PHENOTYPELIST clinvar_NUMSUBMIT clinvar_GUIDELINES)];

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

my ($gnomad_three_file, @extra) = grep { $_->name eq 'custom_gnomadV3_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_gnomadV3_vcf';
}
if($gnomad_three_file) {
    my $vcf_path = $gnomad_three_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'gnomADg';
    $annotation_info->{gnomad_filter} = 'true';
    $annotation_info->{check_existing} = 'true';
    $annotation_info->{vcf_fields} = ['AF','AF_ami','AF_oth','AF_afr','AF_sas','AF_asj','AF_fin','AF_amr','AF_nfe','AF_eas'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}
$inputs->{vep_custom_annotations} = \@vep_custom_annotations;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, $inputs);
