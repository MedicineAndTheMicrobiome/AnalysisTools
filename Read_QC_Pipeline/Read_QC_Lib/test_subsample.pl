#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::Subsample;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $bin=$cfg->val("subsample", "bin");
my $sample_size=10;
#my $sample_size=$cfg->val("subsample", "sample_size");

my $test_data="test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $mod= new Read_QC_Lib::Subsample;

$mod->set_executable_path($bin);
$mod->set_input_fastq($test_data);
$mod->set_sample_size($sample_size);
$mod->set_output_directory("tmp");

$mod->perform_qc();

my $rem_fn=$mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


