#!/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::AdapterTrim;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $adptrim_bin=$cfg->val("adapter_trim", "bin");

my $test_data="../test_data/TAGGCATG-CTCTCTAT_R1.fastq";

my $adptrim_mod= new Read_QC_Lib::AdapterTrim;

$adptrim_mod->set_executable_path($adptrim_bin);
$adptrim_mod->add_adapter_sequence("TAGGCATG");
$adptrim_mod->add_adapter_sequence("CTCTCTAT");
$adptrim_mod->set_temporary_directory($global_scratch);
$adptrim_mod->set_input_fastq($test_data);
$adptrim_mod->set_output_directory("/usr/local/devel/DAS/users/kli/SVN/DAS/Read_QC_Pipeline/Read_QC_Lib/tmp");

$adptrim_mod->perform_qc();

my $rem_fn=$adptrim_mod->get_processed_fastq();
print STDERR "Processed fastq: $rem_fn\n";


