#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::QualityTrim;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $qvtrim_bin=$cfg->val("quality_trim", "bin");

my $qthres=$cfg->val("quality_trim", "quality_threshold");
my $mincut=$cfg->val("quality_trim", "minimum_length");

my $test_data="test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $qvtrim_mod= new Read_QC_Lib::QualityTrim;

$qvtrim_mod->set_executable_path($qvtrim_bin);
$qvtrim_mod->set_quality_threshold($qthres);
$qvtrim_mod->set_minimum_length($mincut);
$qvtrim_mod->set_temporary_directory($global_scratch);
$qvtrim_mod->set_input_fastq($test_data);
$qvtrim_mod->set_output_directory("tmp");

$qvtrim_mod->perform_qc();

my $rem_fn=$qvtrim_mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


