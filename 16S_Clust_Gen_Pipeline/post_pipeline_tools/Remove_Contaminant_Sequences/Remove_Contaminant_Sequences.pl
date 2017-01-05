#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin ();
use Getopt::Std;
use FileHandle;
use File::Basename;
use vars qw($opt_a $opt_g $opt_n $opt_f $opt_l $opt_c $opt_o);

my $PIPELINE_UTIL_PATH="$FindBin::Bin/../../pipeline_utilities";
my $EXTRACT_FASTA_PATH="$PIPELINE_UTIL_PATH/Extract_Record_From_FASTA_By_List.pl";

my $DEF_ID_LEVEL=0.10;

print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";

getopts("a:g:n:f:l:c:o:");

my $usage = "usage: 

$0 
	-a <an.list file>
	-g <groups file (sequence_id / sample_id) >
	-n <names file (sequence_id / seq_id,seq_id,...) >
	-c <list of contaminant/control sample ids>
	-f <fasta file>
	-o <output_file_root>

	[-l <shared identity level, default=$DEF_ID_LEVEL>]

	This script will take in key files from a completed mothur run
	and a list of samples that are considered negative controls, and 
	generate new inputs for another mothur run after removing sequences
	that cluster with sequences from the negative controls.

	Mothur Run Directory:
	or
	<project_id>.unique.good.filter.unique.precluster.pick.an.list
	<project_id>.unique.good.filter.unique.precluster.pick.names

	Format:

	label	numOtus	Otu00001	Otu00002	Otu00003	...	Otu0000N
	unique	#####	<read_ids>	<read_ids>	<read_ids>	...	<read_ids>
	0.01	#####	<read_ids>	<read_ids>	<read_ids>	...	<read_ids>
	0.03	####	<read_ids>	<read_ids>	<read_ids>	...	<read_ids>
	...	
	0.20	###	<read_ids>	<read_ids>	<read_ids>	...	<read_ids>

	Where <read_ids> are a list of comma separated read ids from the fasta file.


	Steps:

	Load Data/Mapping Files:
		1.) Read in the *.an.list file, extracting out row containing the target 'label' matching
			the specified 'shared identity level'.
		2.) Load read -> Representative read mapping
		3.) Load list of contaminant samples
		4.) Load sample to read mapping file (groups file)

	Build OTU Filter:
		5.) For each contaminant sample ID:
			For each read
				Look up representive read
				Look up OTU, and mark as contaminated OTU
		
	Perform Filtering and new Groups file:
		6.) For each read (in names file) 
			Look up represetative read
			Look up OTU, if not marked as contaminated, output reads

	Extract FASTA
		7.) Extract new FASTA file with uncomtaminant associated reads


";

if(!(
	defined($opt_a) && 
	defined($opt_g) && 
	defined($opt_n) && 
	defined($opt_f) && 
	defined($opt_c) && 
	defined($opt_o) 
)){
	die $usage;
}

my $AnListFilename=$opt_a;
my $GroupsFilename=$opt_g;
my $NamesFilename=$opt_n;
my $ContaminantsFilename=$opt_c;
my $FastaFilename=$opt_f;
my $OutputFilenameRoot=$opt_o;
my $LevelofIdentity=$opt_l;

if(!defined($opt_l)){
	$LevelofIdentity=$DEF_ID_LEVEL;
}

$OutputFilenameRoot="$OutputFilenameRoot\." . sprintf("%3.2f",$LevelofIdentity); 


print STDERR "\n";
print STDERR "An.list Filname: $AnListFilename\n";
print STDERR "Groups Filename: $GroupsFilename\n";
print STDERR "Names Filename: $NamesFilename\n";
print STDERR "Contaminant List Filename: $ContaminantsFilename\n";
print STDERR "FASTA Filename: $FastaFilename\n";
print STDERR "Level of Co-mingling: $LevelofIdentity\n";
print STDERR "Output Filename Root: $OutputFilenameRoot\n";
print STDERR "\n";


###############################################################################

sub load_anlist_file{
	my $fname=shift;
	my $cutoff=shift;

	# First line is Otu ID list
	# Example:

	# label   numOtus	<OTU ID1>	<OTU ID2> ... <OTU IDn>
	# unique  80810		<read ids,,,>	<read ids,,,> ... <
	# 0.01    41150		<read ids,,,>	<read ids,,,> ...
	# 0.02    17537		<read ids,,,>	<read ids,,,> ...
	# 0.03    8747		<read ids,,,>	<read ids,,,> ...
	#

	# Read first line
	# Read lines until target cutoff found
	# Return hash

	print STDERR "Loading an.list file: $fname\n";
	print STDERR "Looking for: $cutoff\n";
	open(FH, "<$fname") || die "Could not open $fname.\n";

	# Extract the header row
	my $header=<FH>;
	my @otu_ids=split "\t", $header;

	# Find row with the correct cutoff, then extract those columns.

	my %otu_read_hash; # Mapping from read to OTU ID
	my @avail_levels;
	my $found_level=0;
	my $numOtusLoaded;

	while(<FH>){
		chomp;
		my @fields=split "\t", $_;

		my $label=shift @fields;
		my $numOtus=shift @fields;

		push @avail_levels, $label;
		
		if(($label eq $cutoff) || ($label==$cutoff)){
			print STDERR "Found cutoff: $label...\n";
			print STDERR "Loading OTU Mapping into Memory.\n";
			print STDERR "  Number of OTUs claimed: $numOtus\n";
			my $num_otus_found=$#fields+1;
			print STDERR "  Number of OTUs found: $num_otus_found\n";
			if($num_otus_found!=$numOtus){
				print STDERR "Error: Number of OTUs found ($num_otus_found) not equal to claimed ($numOtus)\n";
				die;
			}
			$numOtusLoaded=$numOtus;

			my $num_represented_reads=0;
			for(my $otu_id=0; $otu_id<$numOtus; $otu_id++){
				my @read_ids_arr=split ",", (shift @fields);
				foreach my $read_id(@read_ids_arr){
					$otu_read_hash{$read_id}=$otu_id;
					$num_represented_reads++;
				}
			}
			print STDERR "Num Represented Reads Found: $num_represented_reads\n";
			$found_level=1;
			last;
		}
	}

	if(!$found_level){
		print STDERR "Available Levels: ", (join ",", @avail_levels), "\n";
		die "Error: Could not find targeted level: $cutoff\n";
	}

	print STDERR "ok.\n\n";

	# Read ID -> OTU ID
	return(\%otu_read_hash, $numOtusLoaded);

}

#------------------------------------------------------------------------------
	
sub load_names_file{
	my $fname=shift;

	# Example:
	# <read_id>	<read ids,,,>
	# <read_id>	<read ids,,,>
	# ...
	# <read_id>	<read ids,,,>

	print STDERR "Loading names file: $fname\n";
        open(FH, "<$fname") || die "Could not open $fname.\n";
	
	my %read_to_rep_hash; # Mapping from read to representative read
	my $tot_reads=0;
	my $tot_represented_reads=0;

	while(<FH>){
		chomp;
		my ($reptv, $reptd)=split "\t", $_;
		my @reptd_arr=split ",", $reptd;
		foreach my $reptd_id(@reptd_arr){
			$read_to_rep_hash{$reptd_id}=$reptv;
			$tot_represented_reads++;
		}
		$tot_reads++;
	}
	
	close(FH);

	print STDERR "Total Representative Reads: $tot_reads\n";
	print STDERR "Total Represented Reads: $tot_represented_reads\n";

	print STDERR "ok.\n\n";

	# Read ID -> Representative Read ID
	return(\%read_to_rep_hash);

}

#------------------------------------------------------------------------------

sub load_groups_file{
	my $fname=shift;

	# Example:
	# <read_id>	<group_id>
	# <read_id>	<group_id>
	# <read_id>	<group_id>
	# ...
	# <read_id>	<group_id>
	
	print STDERR "Loading groups file: $fname\n";
	open(FH, "<$fname") || die "Could not open $fname.\n";
		
	my %groups_hash;
	my $num_read_ids=0;
	while(<FH>){
		chomp;
		my ($read_id, $group_id)=split "\t", $_;
		push @{$groups_hash{$group_id}}, $read_id;
		$num_read_ids++;
	}
	
	close(FH);
	
	print STDERR "Num Groups (Samples): ", (scalar keys %groups_hash), "\n";
	print STDERR "Num Reads: ", $num_read_ids, "\n";

	print STDERR "ok.\n\n";

	# Group ID -> Read IDs Arr
	return(\%groups_hash);

}

#------------------------------------------------------------------------------

sub load_list_file{
	my $fname=shift;
	
	print STDERR "Loading list: $fname\n";
	open(FH, "<$fname") || die "Could not open $fname.\n";
	
	my %list_hash;
	my $num_items=0;
	while(<FH>){
		chomp;
		$list_hash{$_}=1;
		$num_items++;
	}

	close(FH);

	print STDERR "Items loaded: $num_items\n";
	print STDERR "ok.\n\n";
	return(\%list_hash);
}

###############################################################################

# Load read/sample/data files

my $repread_to_otu_hash_ref;
my $read_to_repread_hash_ref;
my $sample_to_read_hash_ref;
my $contam_hash_ref;
my $num_otus;

# Load all the mapping files
($repread_to_otu_hash_ref, $num_otus)=load_anlist_file($AnListFilename, $LevelofIdentity);
$read_to_repread_hash_ref=load_names_file($NamesFilename);
$sample_to_read_hash_ref=load_groups_file($GroupsFilename);
$contam_hash_ref=load_list_file($ContaminantsFilename);

my %contaminated_otu;
my $num_otus_contam;
my $num_contam_samples=scalar keys %{$contam_hash_ref};

my $rarefact_fname="$OutputFilenameRoot\.cont_raref\.tsv";
open(RAREFACTION_FH, ">$rarefact_fname") || die "Could not open $rarefact_fname\n";

print STDERR "Building list of contaminated OTUs...\n";

print STDERR sprintf("%30s %10s %35s\n", "Sample ID", "Num Reads", "Cumulative Contam/Total OTUs");
print RAREFACTION_FH (join "\t", ("Sample ID", "Num Reads in Sample", "Cumulative Contam/Total OTUs")) ."\n";
foreach my $cont_samp_id(sort keys %{$contam_hash_ref}){
	my @cont_reads_arr=@{${$sample_to_read_hash_ref}{$cont_samp_id}};
	my $num_reads=($#cont_reads_arr+1);
	
	foreach my $cont_read(@cont_reads_arr){
		my $rep_read=${$read_to_repread_hash_ref}{$cont_read};
		$contaminated_otu{${$repread_to_otu_hash_ref}{$rep_read}}=1;
	}

	my @contam_otus=keys %contaminated_otu;
	$num_otus_contam=$#contam_otus +1;

	print STDERR sprintf("%30s %10i %35s\n", $cont_samp_id, $num_reads, "$num_otus_contam/$num_otus");
	print RAREFACTION_FH (join "\t", ($cont_samp_id, $num_reads, "$num_otus_contam/$num_otus")) ."\n";
}

close(RAREFACTION_FH);

print STDERR "Percentage of OTUs considered contaminated: ", sprintf("%2.3f",100*$num_otus_contam/$num_otus), "%\n";

print STDERR "ok.\n\n";

#------------------------------------------------------------------------------

print STDERR "Filtering reads in contaminated OTUs...\n";

my @kept_read_groups_arr;
my @kept_reads_arr;

my $sample_breakdown_fh="$OutputFilenameRoot\.sample_brkdwn\.tsv";
open(BRKDWN_FH, ">$sample_breakdown_fh") || die "Could not open $sample_breakdown_fh\n";

print STDERR sprintf("%30s %15s %17s %17s %17s %17s", "[Sample ID]", "[Num Seqs Kept]", "[Num Seqs Total]", "[Prop Seq Kep]", "[Num OTUs Kept]", "[Num OTUs Total]"), "\n";
print BRKDWN_FH (join "\t", ("Sample ID", "Num Seq Kept", "Num Seq Total", "Perc Seq Kept", "Num OTUs Kept", "Num OTUs Total", "Perc OTUs Kept")) . "\n";

my $total_seqs=0;
my $num_total_samples=scalar keys %{$sample_to_read_hash_ref};
my $num_total_contam_seqs=0;
foreach my $samp_id (keys %{$sample_to_read_hash_ref}){
	#if(${$contam_hash_ref}{$samp_id}){
	#	print STDERR "Skipping $samp_id\n";
	#	next;
	#}

	my @reads_in_sample=@{${$sample_to_read_hash_ref}{$samp_id}};
	my $num_total_in_sample=$#reads_in_sample+1;
	$total_seqs+=$num_total_in_sample;
	my $num_kept_in_sample=0;
	my $rejected=0;
	my %kept_otus;
	my %rejc_otus;
	
	foreach my $read_id (@reads_in_sample){
		my $rep_read=${$read_to_repread_hash_ref}{$read_id};
		my $otu=${$repread_to_otu_hash_ref}{$rep_read};
		if(!defined($contaminated_otu{$otu})){
			push @kept_read_groups_arr, "$read_id\t$samp_id";
			push @kept_reads_arr, "$read_id\t$samp_id";
			$num_kept_in_sample++;
			$kept_otus{$otu}=1;
		}else{
			$rejc_otus{$otu}=1;
			$num_total_contam_seqs++;
		}
	}

	my $num_kept_otus=scalar keys %kept_otus;
	my $num_rejc_otus=scalar keys %rejc_otus;
	my $num_tot_otus=$num_kept_otus+$num_rejc_otus;

	print STDERR sprintf("%30s %15i %17i %17.3f %17i %17i", $samp_id, $num_kept_in_sample, $num_total_in_sample, $num_kept_in_sample/$num_total_in_sample, $num_kept_otus, $num_tot_otus), "\n";
	print BRKDWN_FH (join "\t", ($samp_id, $num_kept_in_sample, $num_total_in_sample, sprintf("%3.2f", 100*$num_kept_in_sample/$num_total_in_sample), 
		$num_kept_otus, $num_tot_otus,  sprintf("%3.2f", 100*$num_kept_otus/$num_tot_otus))) ."\n";
}

close(BRKDWN_FH);

#------------------------------------------------------------------------------

my $summary_fn="$OutputFilenameRoot\.stats.tsv";
open(SUMMARY_FH, ">$summary_fn") || die "Could not open $summary_fn.\n";
print SUMMARY_FH "# ContamFile\tNumCtrlSamp\tNumNonCtrlSamp\tClustLevel\tNumContamOTUs\tNumTotOTUs\tNumContamSeq\tNumTotlSeq\tPercOTUsContam\tPercSeqsContam\tNumOTUsKept\tNumSeqsKept\n";
print SUMMARY_FH (join "\t", (
	$ContaminantsFilename,
	$num_contam_samples,
	$num_total_samples-$num_contam_samples,

	$LevelofIdentity,
	$num_otus_contam,
	$num_otus,

	$num_total_contam_seqs,
	$total_seqs,

	sprintf("%3.2f", 100*($num_otus_contam/$num_otus)),
	sprintf("%3.2f", 100*($num_total_contam_seqs/$total_seqs)),

	$num_otus-$num_otus_contam,
	$total_seqs-$num_total_contam_seqs
)) . "\n";
close(SUMMARY_FH);


#------------------------------------------------------------------------------

print STDERR "Outputing new Group and FASTA file.\n";

my $no_contam_groups_fn="$OutputFilenameRoot\.no_ctm.groups";
open(FH_GR, ">$no_contam_groups_fn") || die "Could not open $no_contam_groups_fn for writing.\n";
foreach my $group_pair(@kept_read_groups_arr){
	print FH_GR "$group_pair\n";
}
close(FH_GR);

#------------------------------------------------------------------------------

my $no_contam_reads_fn="$OutputFilenameRoot\.no_ctm.groups";
open(FH_LI, ">$no_contam_reads_fn") || die "Could not open $no_contam_reads_fn for writing.\n";
foreach my $read(@kept_reads_arr){
	print FH_LI "$read\n";
}
close(FH_LI);

my $exec_str="$EXTRACT_FASTA_PATH";


###############################################################################

print STDERR "done.\n";

