#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../Annotation_Lib";
use Getopt::Std;
use vars qw($opt_l $opt_P $opt_A $opt_S $opt_p);
use File::Basename;
use POSIX; # for ceil function

use Config::IniFiles;

getopts("l:P:A:S:p:");
my $usage = "usage: 

$0 

	-l <List of sample ids to process, should be from alignment step>
	[-P <project name, default=sample id fname w/o extension>]

	-A <input annotation result directory, eg. AnnotResults>
	-S <output summary tables directory, eg. SummaryTables>

	-p <Annotation parameter file (mostly database and bin directories)>

	This script should be run after the Run_Read_Annotation.pl pipeline script
	has successfully completed on all the samples.  This script will accumulate
	all the processed data and generate summary tables for taxonomy and function
	at all the cutoffs specified.
	
	List Format (-l option):
		<library/sample name> \\t <library/sample fasta file>

		Each column is tab separated.


";

if(!defined($opt_l)||!defined($opt_A)||!defined($opt_S)||!defined($opt_p)){
	die $usage;
}

# required
my $sample_list=$opt_l;
my $annot_dir=$opt_A;
my $sumtab_dir=$opt_S;
my $parameter_fname=$opt_p;

# optional
my $project_name=$opt_P;

if(!defined($project_name) || $project_name eq ""){
	my ($name, $path, $suffix)=fileparse($sample_list, (".txt", ".tsv", ".lst"));
	$project_name=$name;
}

###############################################################################

print STDERR "\n";
print STDERR "Sample List: $sample_list\n";
print STDERR "Project Name: $project_name\n";
print STDERR "Annotation Dir: $annot_dir\n";
print STDERR "Summary Table Dir: $sumtab_dir\n";
print STDERR "Parameter/Ini File: $parameter_fname\n";
print STDERR "\n";

###############################################################################

sub load_file_list{
	my $list=shift;

	system("dos2unix $list");

	open(IN_FH, "<$list") || die "Could not opne $list\n";
	my @load_file;
	while(<IN_FH>){
		chomp;
		if(substr($_, 0, 1) eq "#"){ next;} # Skip comments
		push @load_file, $_;
	}
	close(IN_FH);	
	return(\@load_file);
}

###############################################################################

sub make_dir{
	my $dir=shift;
	print STDERR "Making $dir.\n";
	if(!(-e $dir)){
		mkdir $dir;
	}else{
		print STDERR "  but $dir exists.\n";
	}
	if(!(-e $dir)){
		die "Could not make or find $dir\n";
	}
}

###############################################################################

sub log_timing_stats{
	my $time_rec_begin_ref=shift;
	my $time_rec_end_ref=shift;
	my $time_output_fn=shift;
	my $name=shift;

	my $start=0;
	my $end=0;
	
	for(my $i=0; $i<4; $i++){
		$start+=${$time_rec_begin_ref}[$i];
		$end+=${$time_rec_end_ref}[$i];
	}

	my $tot_time=$end-$start;

	open(FH, ">>$time_output_fn") || die "Could not open $time_output_fn for appending.\n";

	print FH "$name\t$tot_time\n";	

	close(FH);
}

###############################################################################

my $cfg=Config::IniFiles->new(-file => $parameter_fname);

###############################################################################

sub lookup{
	my $group=shift;
	my $field=shift;

	my $val=$cfg->val($group, $field);
	return($val);
}

sub execute{
        my $desc=shift;
        my $cmd=shift;

        print STDERR "\n";
        print STDERR "*************************************************************************\n";
        print STDERR "*  $desc\n";
        print STDERR "*************************************************************************\n";
        print STDERR "\n";

        $cmd=~s/\s+/ /g;
        my $err=system($cmd);
	if($err!=0){ die "Error in $cmd\n";}
}


###############################################################################

#[annotation_data]
my $ec_mapping=lookup("annotation_data","ec_mapping");
my $go_mapping=lookup("annotation_data","go_mapping");
my $pfam_mapping=lookup("annotation_data","pfam_mapping");
my $tigrfam_mapping=lookup("annotation_data","tigrfam_mapping");

#[configuration]
my $cutoffs=lookup("configuration","cutoffs");
$cutoffs=~s/\s+//;
my @cutoff_arr=split ",", $cutoffs;
my $temporary_directory=lookup("configuration", "temporary_directory");

#[annotation_bin]
my $taxa_levels_to_st_bin=lookup("annotation_bin","taxa_levels_to_st_bin");
my $annot_to_st_bin=lookup("annotation_bin","annot_to_st_bin");
my $rename_st_ids_bin=lookup("annotation_bin","rename_st_ids_bin");

###############################################################################
# assign function type to ID to name map files
my @function_types=("ECIDs", "GOFuncIDs", "GOProcIDs", "TIGRFamIDs", "PFamIDs");
my %function_map;
$function_map{"ECIDs"}=$ec_mapping;
$function_map{"GOFuncIDs"}=$go_mapping;
$function_map{"GOProcIDs"}=$go_mapping;
$function_map{"TIGRFamIDs"}=$pfam_mapping;
$function_map{"PFamIDs"}=$tigrfam_mapping;

my $TAXONOMY_DIR="Taxonomic";
my $FUNCTIONAL_DIR="Functional";
my $TAX_SUM_TAB_LIST_FN="taxa_sample_list.tsv";
my $FUNCT_TARGETS_SUM_TAB_LIST_FN="funct_sample_list.tsv";

make_dir("$sumtab_dir");
make_dir("$sumtab_dir/$TAXONOMY_DIR");
make_dir("$sumtab_dir/$FUNCTIONAL_DIR");

print STDERR "Loading Sample ID list...\n";
my $sample_info_arr_ref=load_file_list($sample_list);
my @sample_ids=@{$sample_info_arr_ref};


###############################################################################
# Build Taxa Summary Tables

# Build list of hghr.taxa_name.tsv files
print STDERR "Building taxa target list...\n";
open(FH, ">$sumtab_dir/$TAXONOMY_DIR/$TAX_SUM_TAB_LIST_FN") || die "Could not open $sumtab_dir/$TAXONOMY_DIR/$TAX_SUM_TAB_LIST_FN\n";
foreach my $samp_id (@sample_ids){
	my $hgh_taxa_fn="$annot_dir/$samp_id/blst.cpsid.trmbl.taxa_est.hghr.taxa_names.tsv";
	print FH "$samp_id\t$hgh_taxa_fn\n";
}
close(FH);
 
#perl ~/git/AnalysisTools/Annotation/TaxaLevels_to_SummaryTable/TaxaLevels_to_SummaryTable.pl \
execute(
"Accumulating taxa info into Summary Tables",
"$taxa_levels_to_st_bin
	-l $sumtab_dir/$TAXONOMY_DIR/$TAX_SUM_TAB_LIST_FN
	-o $sumtab_dir/$TAXONOMY_DIR/$project_name
");

###############################################################################
# Build Annotation Summary Tables

foreach my $cutoff (@cutoff_arr){

	print STDERR "Working on cutoff: $cutoff %\n";
	make_dir("$sumtab_dir/$FUNCTIONAL_DIR/$cutoff");

	# Build target annotation file list
	print STDERR "Building target annotation file list...\n";
	open(FH, ">$sumtab_dir/$FUNCTIONAL_DIR/$cutoff/$FUNCT_TARGETS_SUM_TAB_LIST_FN") || die "Could not open $sumtab_dir/$FUNCTIONAL_DIR/$FUNCT_TARGETS_SUM_TAB_LIST_FN\n";
	foreach my $samp_id (@sample_ids){
		my $hgh_taxa_fn="$annot_dir/$samp_id/$cutoff/blst.cpsid.trmbl.accm.$cutoff.taxa_filt.kept";
		print FH "$samp_id\t$hgh_taxa_fn\n";
	}
	close(FH);

	# Build functional summary tables
	#perl ~/git/AnalysisTools/Annotation/Annotation_to_SummaryTable/Annotation_to_SummaryTable.pl \
	execute(
	"Accumulating annotation into Summary Tables",
	"$annot_to_st_bin
		-l $sumtab_dir/$FUNCTIONAL_DIR/$cutoff/$FUNCT_TARGETS_SUM_TAB_LIST_FN
		-o $sumtab_dir/$FUNCTIONAL_DIR/$cutoff/$project_name\.$cutoff
	");

	print STDERR "Translating/Renaming ID's to display names...\n";
	foreach my $funct_type(@function_types){

		my $map_file=$function_map{$funct_type};
		print STDERR "Translating: $funct_type ($map_file)\n";

		#perl ~/git/AnalysisTools/Profile/SummaryTableUtilities/Rename_Summary_Table_Categories.pl \
		execute(
		"Rename summary table categories from IDs to display names.",
		"$rename_st_ids_bin
			-i $sumtab_dir/$FUNCTIONAL_DIR/$cutoff/$project_name\.$cutoff\.$funct_type\.summary_table.tsv
			-o $sumtab_dir/$FUNCTIONAL_DIR/$cutoff/$project_name\.$cutoff\.$funct_type\.rnmd.summary_table.tsv
			-m $map_file
			-l \";\"
		");
	}

}

###############################################################################

print STDERR "Done.\n";

