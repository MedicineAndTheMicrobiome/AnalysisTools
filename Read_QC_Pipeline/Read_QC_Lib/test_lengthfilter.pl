#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::LengthFilter;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $bin=$cfg->val("length_filt", "bin");

my $min_len=$cfg->val("length_filt", "min_length");
my $max_len=$cfg->val("length_filt", "max_length");

$min_len=230;
$max_len=260;

my $test_data="test_data/TAGGCATG-CTCTCTAT_R1_shortlong.fasta";

my $mod= new Read_QC_Lib::LengthFilter;

$mod->set_executable_path($bin);

$mod->set_min_length($min_len);
$mod->set_max_length($max_len);

$mod->set_temporary_directory($global_scratch);
$mod->set_input_fasta($test_data);
$mod->set_output_directory("tmp");

$mod->perform_qc();

my $rem_fn=$mod->get_processed_fasta();
print STDERR "Processed fastq: $rem_fn\n";


