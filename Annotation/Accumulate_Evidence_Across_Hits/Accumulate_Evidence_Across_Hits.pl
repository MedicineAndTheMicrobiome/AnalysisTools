#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_a $opt_o $opt_h $opt_c);

getopts("a:o:hc:");

my $usage = "
	usage:
	$0
	-a <uniref alignments info (query_id, perc_comp, uniref_id, ... ) >
	-o <output file>
	-c <cutoffs, eg. \"45,60,75,90\">
	[-h (suppress header flag)]

	This script will read in a uniref alignment info file and generate a
	best guess of the annotation for each read at various cutoffs.  The
	lower the cutoff, the more of the low identity annotation will be
	included when choosing the best non-null feature. 

	Input columns:

	1.) Read ID
	2.) Percent Composite Identity
	3.) Read Length
	4.) UniRef ID
	5.) Description
	6.) Taxa ID
	7.) Subject Length
	8.) PFam IDs
	9.) TIGRFam IDs
	10.) GO Process
	11.) GO Function
	12.) GO Component
	13.) EC

	[Unused]
	13.) NA (Subject composite idenity)
	14.) NA (Subject length)
	15.) Alignment Identity
	16.) Alignment Bit Score
	17.) Alignment E-val

	Output selects the best non-null results across the annotations:

	1.) Read ID
	2.) Percent Composite Identity
	3.) Length
	4.) Description
	5.) PFam IDs
	6.) TIGRFam IDs
	7.) GO Process
	8.) GO Function
	9.) GO Component 
	10.) EC

";

###############################################################################

if(!defined($opt_a) || !defined($opt_o) || !defined($opt_c)){
	die $usage;
}

my $alignments_file=$opt_a;
my $output_file=$opt_o;
my $suppress_header=($opt_h eq "1");
my @cutoffs=split ",", $opt_c;

###############################################################################

#my @cutoffs=(45, 60, 75, 90);
my $num_cutoffs=$#cutoffs+1;
my %fh_hash;

foreach my $cutoff(@cutoffs){
	my $ext=sprintf("%2.0f", $cutoff);
	$ext=~s/^0\.//;
	my $fname="$output_file\.$ext\.tsv";
	$fh_hash{$cutoff}=FileHandle->new;
	$fh_hash{$cutoff}->open(">$fname");
	if(!$fh_hash{$cutoff}){
		die "Could not open $fname.\n";
	}
}

if(!$suppress_header){
	print STDERR "Header being written.\n";
	my $fields_str;
	$fields_str=join "\t", (
		"ReadID",
		"PercCompID",
		"Length",
		"Description",
		"PFamIDs",
		"TIGRFamIDs",
		"GOProcIDs",
		"GOFuncIDs",
		"GoCompIDs",
		"ECIDs"	
	);

	foreach my $cutoff(@cutoffs){
		print {$fh_hash{$cutoff}} "# $fields_str\n";
	}
}else{
	print STDERR "Header for output file has been suppressed.\n";
}

###############################################################################
# Input column definition

my $ID_COL=0;
my $PERC_ID_COL=1;
my $READ_LEN=2;
my $UNIREF_COL=3;
my $DESC_COL=4;
my $TAXA_COL=5;
my $LEN_COL=6;
my $PFAM_COL=7;
my $TIGRFAM_COL=8;
my $GO_PROC_COL=9;
my $GO_FUNC_COL=10;
my $GO_COMP_COL=11;
my $EC_COL=12;

###############################################################################

sub print_rec{
	my $recs=shift;
	foreach my $rec(@{$recs}){
		my $rec_str=join ",", @{$rec};
		print "$rec_str\n";
	}
}

sub get_perc_ids{
	my $recs=shift;
	my @perc_ids;

	foreach my $rec(@{$recs}){
		push @perc_ids, ${$rec}[$PERC_ID_COL];
        }
	return(\@perc_ids);
}
		
sub get_best_function{
	my $rec_ref=shift;

	my ( $best_perc_id, $best_len, $best_desc, 
		$best_pfam, $best_tigrfam, 
		$best_go_proc, $best_go_func, $best_go_comp, $best_ec)=
		(0,"","","","","","","","");

	# These records are already assumed to be sorted by increasing order
	foreach my $rec(@{$rec_ref}){
		
		my ($perc_id, $len, $desc, $pfam, $tigrfam, $go_proc, $go_func, $go_comp, $ec)=(
			${$rec}[$PERC_ID_COL],
			${$rec}[$LEN_COL],
			${$rec}[$DESC_COL],
			${$rec}[$PFAM_COL],
			${$rec}[$TIGRFAM_COL],
			${$rec}[$GO_PROC_COL],
			${$rec}[$GO_FUNC_COL],
			${$rec}[$GO_COMP_COL],
			${$rec}[$EC_COL]
		);

		# If there is a better hit, overwrite the old hit.
		if($perc_id >= $best_perc_id){
			$best_perc_id=$perc_id;
			if($desc ne ""){

				if($best_desc eq ""){
					$best_desc=$desc;
				}elsif(!($desc=~/Uncharacterized protein/ ||
				   $desc=~/Putative uncharacterized protein/)){
					$best_desc=$desc;	
				}else{
					# Ignore
				}
			}
			if($len ne "NA"){ $best_len=$len; }
			if($pfam ne "NA"){ $best_pfam=$pfam; }
			if($tigrfam ne "NA"){ $best_tigrfam=$tigrfam; }
			if($go_proc ne "NA"){ $best_go_proc=$go_proc; }
			if($go_func ne "NA"){ $best_go_func=$go_func; }
			if($go_comp ne "NA"){ $best_go_comp=$go_comp; }
			if($ec ne "NA"){ $best_ec=$ec; }
		}
	}	

	return(($best_perc_id, $best_len, $best_desc,
                $best_pfam, $best_tigrfam,
                $best_go_proc, $best_go_func, $best_go_comp, $best_ec));
	
}

sub randomize{
	my $arr_ref=shift;
	my @rnd_idx;
	my @rnd_val;
	my $num_rec=$#{$arr_ref}+1;
	
	for(my $i=0; $i<$num_rec; $i++){
		push @rnd_idx, rand;
	}

	my @sort_ix=sort {$rnd_idx[$a] <=> $rnd_idx[$b]} 0..$#{$arr_ref};

	for(my $i=0; $i<$num_rec; $i++){
		push @rnd_val, ${$arr_ref}[$sort_ix[$i]];
	}
	
	return(\@rnd_val);
}

sub process_records{

	my $rec_arr_ref=shift;

	if($#{$rec_arr_ref}==-1){
		# Skip if record is empty.
		return;
	}
	
	my $read_id=${${$rec_arr_ref}[0]}[$ID_COL];

	#print "Before RND:\n";
	#print_rec($rec_arr_ref);

	# Randomize so ties will be in random order
	$rec_arr_ref=randomize($rec_arr_ref);

	#print "After RND:\n";
	#print_rec($rec_arr_ref);
	
	my @perc_id_ref=@{get_perc_ids($rec_arr_ref)};

	# Sort in ascending order, so best is seen last
	my @sort_ix=sort {$perc_id_ref[$a] <=> $perc_id_ref[$b]} 0..$#perc_id_ref;

	# Apply sort
	my @sorted_records;
	for(my $i=0; $i<=$#sort_ix; $i++){
		push @sorted_records, ${$rec_arr_ref}[$sort_ix[$i]];
	}

	#print "Sorted:\n";
	#print_rec(\@sorted_records);
	#print "------------------------------\n";

	# Generate groups of records based on cutoffs
	my %grouped_hash;
	for(my $i=0; $i<=$#sort_ix; $i++){
		my $rec_ref=$sorted_records[$i];
		my $perc_id=${$rec_ref}[$PERC_ID_COL];

		foreach my $cutoff(@cutoffs){
			if($perc_id >= $cutoff){
				push @{$grouped_hash{$cutoff}}, $rec_ref;
			}
		}
	}

	# Analyse each group
	foreach my $cutoff(@cutoffs){
		#print "Group at $cutoff:\n";
		#print_rec($grouped_hash{$cutoff});
		my @best=get_best_function($grouped_hash{$cutoff});

		if($best[0]!=0){
			my $best_str=join "\t", @best;
			print {$fh_hash{$cutoff}} "$read_id\t$best_str\n";
		}
	}

}

###############################################################################
# Process alignments file


#open(FH, "<$alignments_file") || die "Could not open $alignments_file\n";
open(FH, "sort -k 1,1 -k2,2nr $alignments_file |") || die "Could not open $alignments_file\n";

my @records=();
my $last_rec="";

my $num_incomplete_records=0;
while(<FH>){
	chomp;

	my @fields=split "\t", $_;

	my $num_fields=$#fields+1;
	if($num_fields<12){
		$num_incomplete_records++;
		#print STDERR "Incomplete Record: '$_'\n";
		next;
	}
	
	if(($fields[$ID_COL] ne $last_rec)){
		process_records(\@records);
		@records=();
	}

	push @records, \@fields;

	$last_rec=$fields[$ID_COL];
}
process_records(\@records);

close(FH);

print STDERR "Number of Incomplete Records Discarded: $num_incomplete_records\n";

###############################################################################

print STDERR "Done.\n";

