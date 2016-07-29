#!/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use Read_QC_Lib::Dust;
use Config::IniFiles;


my $cfg=Config::IniFiles->new(-file => "../read_qc.ini"); 
my $global_scratch=$cfg->val("GLOBAL", "scratch_dir");
my $dust_bin=$cfg->val("dust_screen", "bin");
my $cutoff=$cfg->val("dust_screen", "cutoff");

my $test_data="../test_data/TAGGCATG-CTCTCTAT_R1.fasta";

my $dust_mod= new Read_QC_Lib::Dust;

$dust_mod->set_executable_path($dust_bin);
$dust_mod->set_cutoff($cutoff);
$dust_mod->set_temporary_directory($global_scratch);
$dust_mod->set_input_fasta($test_data);
$dust_mod->set_output_directory("/usr/local/devel/DAS/users/kli/SVN/DAS/Read_QC_Pipeline/Read_QC_Lib/tmp");

$dust_mod->perform_qc();

my $rem_fn=$dust_mod->get_remove_list_path();
print STDERR "Remove List: $rem_fn\n";


