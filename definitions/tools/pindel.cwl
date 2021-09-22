#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
arguments: ["/usr/bin/perl", "pindel_helper.pl"]
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.4.2"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'pindel_helper.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use IO::File;

            my $normal_bam = pop @ARGV;
            if (defined $normal_bam) {
                unless (@ARGV > 5) {
                    die "Usage: $0 normal.bam tumor.bam insert_size normal_sample_name tumor_sample_name <args>";
                }
                my ($tumor_bam, $insert_size, $normal_name, $tumor_name, @args) = @ARGV;
                my $fh = IO::File->new("> pindel.config");

                $fh->say(join("\t", $normal_bam, $insert_size, $normal_name));
                $fh->say(join("\t", $tumor_bam, $insert_size, $tumor_name));
                $fh->close;

                exit system(qw(/usr/bin/pindel -i pindel.config -w 30 -T 4 -o all), @args);
            } else {
                unless (@ARGV > 4) {
                    die "Usage: $0 tumor.bam insert_size normal_sample_name tumor_sample_name <args>";
                }
                my ($tumor_bam, $insert_size, $normal_name, $tumor_name, @args) = @ARGV;
                my $fh = IO::File->new("> pindel.config");

                $fh->say(join("\t", $tumor_bam, $insert_size, $tumor_name));
                $fh->close;

                exit system(qw(/usr/bin/pindel -i pindel.config -w 30 -T 4 -o all), @args);
            }

inputs:
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
        inputBinding:
            position: 1
    normal_bam:
        type: File?
        secondaryFiles: ["^.bai"]
        inputBinding:
            position: 20
    insert_size:
        type: int
        default: 400
        inputBinding:
            position: 2
    tumor_sample_name:
        type: string
        default: 'TUMOR'
        inputBinding:
            position: 4
    normal_sample_name:
        type: string
        default: 'NORMAL'
        inputBinding:
            position: 3
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-f"
            position: 5
    chromosome:
        type: string?
        inputBinding:
            prefix: "-c"
            position: 6
    region_file:
        type: File?
        inputBinding:
            prefix: "-j"
            position: 7
outputs:
    deletions:
        type: File
        outputBinding:
            glob: "all_D"
    insertions:
        type: File
        outputBinding:
            glob: "all_SI"
    tandems:
        type: File
        outputBinding:
            glob: "all_TD"
    long_insertions:
        type: File
        outputBinding:
            glob: "all_LI"
    inversions:
        type: File
        outputBinding:
            glob: "all_INV"
