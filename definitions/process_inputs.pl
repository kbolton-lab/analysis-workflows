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
			reference => { input_type => 'File' },
			UMI_paired => { input_type => 'Number' },
			max_read_error_rate => { input_type => 'Number' },
			max_base_error_rate => { input_type => 'Number' },
			min_base_quality => { input_type => 'Number' },
			max_no_call_fraction => { input_type => 'Number' },
		],
	};
}

my @inputs = $build->inputs;
my $input_processor = InputProcessor->get($build->id);
my $inputs = $input_processor->simple_inputs;

# Store Sample_Name (Key) with the Build Sample Name
$inputs->{sample_name} = $build->subject->name;

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

my @min_reads_input = grep { $_->name eq 'min_reads' } @inputs;
my @min_reads = map { 0 + $_->value_id } @min_reads_input;

$inputs->{min_reads} = \@min_reads;

my @read_structure_inputs = Genome::Model::Build::Input->get(build_id => $build->id, 'name LIKE' => 'read_structure_R%', -order_by => 'name');
my @read_structure = map $_->value_id, @read_structure_inputs;

$inputs->{read_structure} = \@read_structure;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, $inputs);
