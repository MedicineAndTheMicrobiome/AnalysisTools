#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::AbsoluteTrim;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $trim_bin=$cfg->val("absolute_trim", "bin");

my $trim_start=$cfg->val("absolute_trim", "keep_start_1based");
my $trim_end=$cfg->val("absolute_trim", "keep_end_1based");

my $test_data="test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $trim_mod= new Read_QC_Lib::AbsoluteTrim;

$trim_mod->set_executable_path($trim_bin);

$trim_mod->set_trim_start_keep($trim_start);
$trim_mod->set_trim_end_keep($trim_end);

$trim_mod->set_temporary_directory($global_scratch);
$trim_mod->set_input_fastq($test_data);
$trim_mod->set_output_directory("tmp");

$trim_mod->perform_qc();

my $rem_fn=$trim_mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


