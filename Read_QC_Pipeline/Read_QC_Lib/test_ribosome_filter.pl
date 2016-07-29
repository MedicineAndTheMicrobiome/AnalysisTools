#!/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::RibosomeFilter;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $bin=$cfg->val("ribosomal_screen", "bin");

my $test_data="../test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $mod= new Read_QC_Lib::RibosomeFilter;

$mod->set_executable_path($bin);

$mod->set_alignment_identity_threshold($cfg->val("ribosomal_screen", "alignment_identity_threshold"));
$mod->set_alignment_coverage_threshold($cfg->val("ribosomal_screen", "alignment_coverage_threshold"));
$mod->set_alignment_length_threshold($cfg->val("ribosomal_screen", "alignment_length_threshold"));
$mod->set_database_list($cfg->val("ribosomal_screen", "database_list"));

$mod->set_temporary_directory($global_scratch);
$mod->set_input_fastq($test_data);
$mod->set_output_directory("/usr/local/devel/DAS/users/kli/SVN/DAS/Read_QC_Pipeline/Read_QC_Lib/tmp");

$mod->perform_qc();

my $rem_fn=$mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


