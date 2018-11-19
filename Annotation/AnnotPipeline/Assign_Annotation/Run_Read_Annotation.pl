#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../Annotation_Lib";
use Getopt::Std;
use vars qw($opt_l $opt_P $opt_I $opt_O $opt_p $opt_c $opt_s);
use File::Basename;
use POSIX; # for ceil function
use Config::IniFiles;

getopts("l:P:I:O:p:c:s:");
my $usage = "usage: 

$0 

	-l <List of sample ids to process, should be from alignment step>
	[-P <project name, default=sample id fname w/o extension>]
	-I <Input Directory:  read-to-uniprot blastx results, should be output directory of alignment step>
		(Input file will blast.out)
	-O <Result/root directory of blast output, should be output directory of alignment step>
		(output will be saved here, input file will blast.out)

	-p <Annotation parameter file (mostly database and bin directories)>
	[-c <Host/Contaminate taxa filter file, (taxa_id, taxa_name)-tuple>]

	[-s <Parallelization offset, default=1>,<Multiplier, default=1>]

	This script will read in a results directory which contains a composite identity
	file per sample ID and associate each read with annotation and accumulate it
	for each sample.  When completed, the Generate_Summary_Tables.pl script can 
	be applied to accumulate all the samples together.

	
	List Format (-l parameter):
		<library/sample name> \\t <library/sample fasta file>

		Each column is tab separated.


	Input directory (-I parameter):
		The directory will contains a subdirectory for each sample id in the
		list file.  Each sample will have a blast output.

		<input_directory>/<sample_id>/blast.out

		Intermediate output will look like:
		<input_directory>/<sample_id>/blast.out.<steps>

		Output summary tables will be in the <input_directory>:
			<input_directory>/taxonomic_summary_tables
			<input_directory>/functional_summary_tables

	Offset (-s option):
		If the -s option is used, then only a subset of the reads
		will be processed. If offset=1, and multipler=1, then
		items 1,2,3, ... N sample will be processed.  If offset=3 and
		mulitplier=2, then 3,5,7, ... N sample will be processed.
		Offset can be thought of the job ID starting from 1, and multiplier
		is the number of processes to run in parallel.  


	
";

if(!defined($opt_l)||!defined($opt_I)||!defined($opt_O)||!defined($opt_p)){
	die $usage;
}

# required
my $sample_list=$opt_l;
my $input_dir=$opt_I;
my $output_dir=$opt_O;
my $parameter_fname=$opt_p;

# optional
my $project_name=$opt_P;
my $offset_param="1,1"; #opt_s
my $contam_file=""; # opt_c

if(defined($opt_s)){
	$offset_param=$opt_s;	
}
my ($offset, $multiplier)=split ",", $offset_param;

if(defined($opt_c)){
	$contam_file=$opt_c;
}

if(!defined($project_name) || $project_name eq ""){
	$project_name=$sample_list;
	$project_name=~s/\.txt$//;
	$project_name=~s/\.tsv$//;
	$project_name=~s/\.lst$//;
}

###############################################################################

print STDERR "\n";
print STDERR "Sample List: $sample_list\n";
print STDERR "Project Name: $project_name\n";
print STDERR "Input Blast Dir: $input_dir\n";
print STDERR "Output Annotation Dir: $output_dir\n";
print STDERR "Parameter/Ini File: $parameter_fname\n";
print STDERR "Subset Paramters: Offset=$offset Multiplier=$multiplier\n";
print STDERR "Contam/Host Filter: $contam_file\n";
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
		print STDERR "$dir exists.\n";
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
my $sample_info_arr_ref=load_file_list($sample_list);

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
my $TrEMBL_annotation=lookup("annotation_data","TrEMBL_annotation");
my $ncbi_taxa_nodes=lookup("annotation_data","ncbi_taxa_nodes");
my $ncbi_taxa_names=lookup("annotation_data","ncbi_taxa_names");
my $ncbi_taxa_levels=lookup("annotation_data","ncbi_taxa_levels");

#[configuration]
my $cutoffs=lookup("configuration","cutoffs");
$cutoffs=~s/\s+//;
my @cutoff_arr=split ",", $cutoffs;
my $temporary_directory=lookup("configuration", "temporary_directory");

#[annotation_bin]
my $perc_comp_id_bin=lookup("annotation_bin","perc_comp_id_bin");
my $join_aln_w_clst_anno_bin=lookup("annotation_bin","join_aln_w_clst_anno_bin");

my $estimate_taxa_bin=lookup("annotation_bin","estimate_taxa_bin");
my $get_high_taxa_levels_bin=lookup("annotation_bin","get_high_taxa_levels_bin");

my $accumulate_annot_bin=lookup("annotation_bin","accumulate_annot_bin");
my $filter_annot_by_taxa_bin=lookup("annotation_bin","filter_annot_by_taxa_bin");

###############################################################################

my $offset_str_width=ceil(log($multiplier)/log(10));
my $offset_str=sprintf("%0$offset_str_width"."i", $offset);

#my $timing_log_fn="$result_dir/timing_log.$offset_str.tsv";

my $num_records=$#{$sample_info_arr_ref}+1;
$offset--; # Start offset from 1 less than specified, so we start from 0, instead of 1

my $blast_out="blast.out";

make_dir($output_dir);


for(my $idx=$offset; $idx<$num_records; $idx+=$multiplier){

	print STDERR "#########################################################################################\n";
	print STDERR "\nWorking on record: $idx.\n";
	my $samp_info=${$sample_info_arr_ref}[$idx];
	my ($samp_name)=split "\t", $samp_info;
	
	print STDERR "Sample Name: $samp_name\n\n";
	print STDERR "#########################################################################################\n";

	make_dir("$output_dir/$samp_name");
	execute(
	"Linking input blast.out.comp_id results to output directory...",
	"ln -s $input_dir/$samp_name/$blast_out\.comp_id $output_dir/$samp_name/blst.cpsid"
	);

	#----------------------------------------------------------------------

	#perl ~/git/AnalysisTools/Column/remap_column_values.pl \
	execute(
	"Joining composite alignment results with TrEMBL via UniRef ID",
	"$join_aln_w_clst_anno_bin 
		-i $output_dir/$samp_name/blst.cpsid
		-m $TrEMBL_annotation 
		-p 
		-c 3 
		-f > 
		$output_dir/$samp_name/blst.cpsid.trmbl 
	"
	);

	#perl ~/git/AnalysisTools/Annotation/Accumulate_Evidence_Across_Hits/Accumulate_Evidence_Across_Hits.pl \
	execute(
	"Accumulating evidence across all annotation/hits for each read",
	"$accumulate_annot_bin 
		-a $output_dir/$samp_name/blst.cpsid.trmbl 
		-o $output_dir/$samp_name/blst.cpsid.trmbl.accm
		-c \"$cutoffs\"
	");
	# Accumulates at various cutoffs, eg. at 45, 60, 75 and 90.
		
	#perl ~/git/AnalysisTools/Annotation/Estimate_Taxa_from_Alignment/Estimate_Taxa_from_Alignment.pl \
	execute(
	"Estimating greatest common ancestor taxonomy across annotation",
	"$estimate_taxa_bin 
		-a $output_dir/$samp_name/blst.cpsid.trmbl 
		-A \"1,2,6\" 
		-t $ncbi_taxa_nodes 
		-n $ncbi_taxa_names 
		-o $output_dir/$samp_name/blst.cpsid.trmbl.taxa_est
	");

	#perl ~/git/AnalysisTools/Annotation/Estimate_Taxa_from_Alignment/Get_Higher_Levels.pl \
	execute(
	"Looking up higher level taxonomic classifications for each read",
	"$get_high_taxa_levels_bin 
		-t $output_dir/$samp_name/blst.cpsid.trmbl.taxa_est 
		-T \"1,7\" 
		-m $ncbi_taxa_names 
		-d $ncbi_taxa_nodes 
		-l $ncbi_taxa_levels 
		-o $output_dir/$samp_name/blst.cpsid.trmbl.taxa_est.hghr
	");

	foreach my $cutoff (@cutoff_arr){

		print STDERR "Working on cutoff: $cutoff %\n";
		make_dir("$output_dir/$samp_name/$cutoff");

		`mv $output_dir/$samp_name/blst.cpsid.trmbl.accm.$cutoff\.tsv $output_dir/$samp_name/$cutoff`;

		if($contam_file ne ""){

			#perl ~/git/AnalysisTools/Annotation/Filter_Annotation_By_TaxaID/Filter_Annotation_By_TaxaID.pl \
			execute(
			"Filtering reads by assigned taxa",
			"$filter_annot_by_taxa_bin 
				-t $output_dir/$samp_name/blst.cpsid.trmbl.taxa_est.hghr.taxa_ids.tsv 
				-a $output_dir/$samp_name/$cutoff/blst.cpsid.trmbl.accm.$cutoff\.tsv 
				-f $contam_file 
				-o $output_dir/$samp_name/$cutoff/blst.cpsid.trmbl.accm.$cutoff\.taxa_filt
			");
		}

	}

	`touch "$output_dir/$samp_name/completed."`;


}

#my @begin=(0,0,0,0);
#my @time_final_end=times;
#my ($name, $path)=fileparse($read_list);
#log_timing_stats(\@begin, \@time_final_end, $timing_log_fn, "Total Execution\t$name\t");

###############################################################################

print STDERR "Done.\n";

