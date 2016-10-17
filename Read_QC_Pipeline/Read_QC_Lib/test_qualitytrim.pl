#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::QualityTrim;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $qvtrim_bin=$cfg->val("quality_trim", "bin_trim");
my $qvfilt_bin=$cfg->val("quality_trim", "bin_filt");

my $trim_qthres=$cfg->val("quality_trim", "trim_quality_threshold");
my $trim_mincut=$cfg->val("quality_trim", "trim_minimum_length");
my $filt_qthres=$cfg->val("quality_trim", "filt_quality_threshold");
my $filt_minperc=$cfg->val("quality_trim", "filt_percent_above");

my $test_data="test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $qvtrim_mod= new Read_QC_Lib::QualityTrim;

$qvtrim_mod->set_trim_executable_path($qvtrim_bin);
$qvtrim_mod->set_filt_executable_path($qvfilt_bin);

$qvtrim_mod->set_trim_quality_threshold($trim_qthres);
$qvtrim_mod->set_trim_minimum_length($trim_mincut);
$qvtrim_mod->set_filt_quality_threshold($filt_qthres);
$qvtrim_mod->set_filt_percent_above($filt_minperc);

$qvtrim_mod->set_temporary_directory($global_scratch);
$qvtrim_mod->set_input_fastq($test_data);
$qvtrim_mod->set_output_directory("tmp");

$qvtrim_mod->perform_qc();

my $rem_fn=$qvtrim_mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


