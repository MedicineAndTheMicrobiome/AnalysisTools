#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use Sys::Hostname;
use File::Basename;
use FindBin ();
use vars qw($opt_i $opt_r $opt_o $opt_d $opt_k $opt_e $opt_l $opt_m $opt_p $opt_c);

my $BLASTALL_BIN=`which blastall`;
my $REVERSE_COMPLEMENT_FASTA_BIN="$FindBin::Bin/Reverse_Complement_FASTA.pl";
my $EXTRACT_RECORD_FASTA_BIN="$FindBin::Bin/Extract_Record_From_FASTA_By_List_Using_STDIO.pl";

my $DEF_BLAST_EVAL=1e-100;

getopts("i:r:o:d:k:e:l:m:p:c");
my $usage = "usage: 
$0 
	-i <input fasta>
	-r <reference fasta, this should already be formatdb'd>
	-o <output fasta>

	[-d <output reject/delete list>]
	[-k <output keep list>]
	[-e <output orient list>]

	[-l <best HSP filter percentage, eg 95>]
	[-m <e-value for blast filtering, default=$DEF_BLAST_EVAL>]
	[-p <blast program, default blastx>]
	[-c <include reference fasta in output fasta>]

	This program will go through all the input sequences, compare them against the 
	reference, in order to determine if the input sequences are in the correct orientation.

	The -l option can be used to filter out input sequences that are not at least the percentage
	specified consecutively.

	This program depends on:

		$BLASTALL_BIN
	and
		$REVERSE_COMPLEMENT_FASTA_BIN
	and 
		$EXTRACT_RECORD_FASTA_BIN

";

if(!(
	defined($opt_i) ||
	defined($opt_r) ||
	defined($opt_o))){
	die $usage;
}

my $input_fasta=$opt_i;
my $reference_fasta=$opt_r;
my $output_fasta=$opt_o;

my $filter_hsp_length=defined($opt_l);
my $filter_cutoff=$opt_l;

my $keep_rejectlist=$opt_d;
my $keep_keeplist=$opt_k;
my $keep_reorientlist=$opt_e;

my $include_reference=defined($opt_c);

my $blast_program="blastx";
if(defined($opt_p)){
	$blast_program=$opt_p;
}

my $blast_eval=$DEF_BLAST_EVAL;
if(defined($opt_m)){
	$blast_eval=$opt_m;
}

# Echo parameters
print STDERR "Input Fasta: $input_fasta\n";
print STDERR "Output Fasta: $output_fasta\n";
print STDERR "Reference: $reference_fasta\n";
print STDERR "Using '$blast_program' to check against reference.\n";
print STDERR "BLAST e-value setting: $blast_eval\n";
if($filter_hsp_length){
	print STDERR "Filtering HSP lengths < $filter_cutoff%\n";
}

# Get unique temporary file name
my $hostname=hostname();
my $pid=$$;
my ($progname)=fileparse($0);
my $tmp_root="$output_fasta\.$progname\.$hostname\.$pid";

# Set names of intermediate results
my $blast_results_file="$tmp_root\.blastout";
my $reorient_list="$tmp_root\.reorient";
my $keep_list="$tmp_root\.keep_list";
my $post_orient_fasta="$tmp_root\.post_orient";

###############################################################################


# Get reference sequence lengths
my $ref_length_hash;
my $ref_sequence_hash;
($ref_length_hash, $ref_sequence_hash)=get_sequence_info($reference_fasta);

# Run blast
run_blast($input_fasta, $reference_fasta, $blast_program, $blast_results_file);

# Determine reorient and filter list
my ($reorient_hash_ref, $filter_hash_ref, $reject_hash_ref, $used_reference_hash_ref)=
	make_reorient_and_filter_list($blast_results_file, $ref_length_hash, $filter_cutoff);

# Perform reorientation
reorient($input_fasta, $reorient_hash_ref, $output_fasta);

# Perform filtering
if($filter_hsp_length){
	print STDERR `mv $output_fasta $post_orient_fasta`;
	filter($post_orient_fasta, $filter_hash_ref,$output_fasta);
	print STDERR `rm $post_orient_fasta`;
}

# Append reference to output if necessary
if($include_reference){
	print STDERR "Including reference in $output_fasta\n";
	open(OUT_FH, ">>$output_fasta") || die "Could not open $output_fasta\n";
	foreach my $id(sort keys %{$used_reference_hash_ref}){
		print STDERR "Including: $id\n";
		print OUT_FH ">$id\n";
		print OUT_FH "${$ref_sequence_hash}{$id}\n";
	}
	close(OUT_FH);
}

if(defined($keep_keeplist)){
	exec_cmd("mv $keep_list $keep_keeplist");
}
if(defined($keep_reorientlist)){
	exec_cmd("mv $reorient_list $keep_reorientlist");
}
if(defined($keep_rejectlist)){
	open(RJ_FH, ">$keep_rejectlist") || die "Could not open $keep_rejectlist\n";
	foreach my $id(sort keys %{$reject_hash_ref}){
		print RJ_FH "$id\t${$reject_hash_ref}{$id}\n";
	}
	close(RJ_FH);
}



###############################################################################

sub exec_cmd{
	my $cmd=shift;
	$cmd=~s/[\t\n]+/ /g;
	print STDERR "Going to execute: $cmd\n";
	system($cmd);
}

###############################################################################

sub get_sequence_info{
	my $file=shift;
	my %length_hash;	
	my %sequence_hash;

	open(FH, "<$file") || die "Could not open $file for reading.\n";

	my $sequence="";
	my $defline;
	my $prev_defline;
	while(<FH>){
		chomp;

		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, \$sequence, \%length_hash, \%sequence_hash);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	process_record($prev_defline, \$sequence, \%length_hash, \%sequence_hash);

	close(FH);

	return(\%length_hash, \%sequence_hash);
}

#------------------------------------------------------------------------------

sub process_record{
	my $defline=shift;
	my $sequence_ref=shift;
	my $length_hash_ref=shift;
	my $sequence_hash_ref=shift;
	
	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Error parsing deflined: $defline\n";
	}

	${$length_hash_ref}{$id}=length(${$sequence_ref});
	${$sequence_hash_ref}{$id}=${$sequence_ref};
}

###############################################################################

sub run_blast{
	my $input_file=shift;
	my $reference_fasta=shift;
	my $blast_program=shift;
	my $results_file=shift;

	my $cmd="
		$BLASTALL_BIN 
			-p $blast_program 
			-i $input_file 
			-d $reference_fasta
			-m 8 
			-e $blast_eval
			-o $results_file
	";
	exec_cmd($cmd);

}

###############################################################################

sub make_reorient_and_filter_list{
	my $results_file=shift;
	my $ref_len_hash_ref=shift;
	my $filter_cutoff=shift;

	open(FH, "<$results_file") || die "Could not open $results_file for reading.\n";

	my %best_hit_hash;
	my %best_hit_ori_hash;
	my %best_hit_align_len_hash;
	my %best_hit_subject_hash;

	while(<FH>){
		chomp;
		
		# Parse the blast results
		my 
			($query_id,$subject_id,$perc_identity,$alignment_length,
			$mismatches,$gap_openings,$q_start,$q_end,$s_start,$s_end,
			$evalue,$bit_score)=split /\t/, $_;


print STDERR "$_\n";
print STDERR "($query_id,$subject_id,$perc_identity,$alignment_length,
                        $mismatches,$gap_openings,$q_start,$q_end,$s_start,$s_end,
                        $evalue,$bit_score)\n";
print STDERR "$query_id: $q_start,$q_end,$s_start,$s_end\n";
				
		if(!defined($best_hit_hash{$query_id})){
			# If sequence not seen yet, take first as best/only hit
			$best_hit_hash{$query_id}=$bit_score;
			$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
			$best_hit_align_len_hash{$query_id}=abs($s_start-$s_end);
			$best_hit_subject_hash{$query_id}=$subject_id;
		}else{
			# Replace best if we find something with a higher bit score
			my $prior_score=$best_hit_hash{$query_id};
			if($prior_score<$bit_score){
				$best_hit_hash{$query_id}=$bit_score;
				$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
				$best_hit_align_len_hash{$query_id}=abs($s_start-$s_end);
				$best_hit_subject_hash{$query_id}=$subject_id;
			}
		
		}

	}

	close(FH);

	# Grab all the query_id's that we need to reorient
	my %reorient_hash;
	foreach my $query_id(keys %best_hit_hash){
		if($best_hit_ori_hash{$query_id}==0){
			$reorient_hash{$query_id}=1;
		}
	}

	# Grab all the query_id's that are longer than cutoff
	my %filter_keep_hash;
	my %reject_hash;
	if(defined($filter_cutoff)){
		my $decimal_cutoff=$filter_cutoff/100.0;
		foreach my $query_id(keys %best_hit_hash){

			#print STDERR "align: $best_hit_align_len_hash{$query_id}\n";
			#print STDERR "subj: ${$ref_len_hash_ref}{$best_hit_subject_hash{$query_id}}\n";

			my $perc_aligned=$best_hit_align_len_hash{$query_id}/${$ref_len_hash_ref}{$best_hit_subject_hash{$query_id}};
			if($perc_aligned>=$decimal_cutoff){
				$filter_keep_hash{$query_id}=1;	
			}else{
				$reject_hash{$query_id}=$perc_aligned;
			}
		}
	}

	# Build unique list of subjects that were used
	my %used_reference_hash;
	foreach my $query_id(keys %best_hit_hash){
		$used_reference_hash{$best_hit_subject_hash{$query_id}}=1;
	}

	# Return a hash of query id's that are in different orientation with reference/subject
	return(\%reorient_hash, \%filter_keep_hash, \%reject_hash, \%used_reference_hash);
	
}

###############################################################################

sub reorient{
	# Calls outside script to reverse complement based on a list.

	my $input_fasta=shift;
	my $reorient_hash_ref=shift;
	my $output_fasta=shift;

	open(FH, ">$reorient_list") || die "Could not open $reorient_list\n";
	foreach my $id(sort keys %{$reorient_hash_ref}){
		print FH "$id\n";
	}
	close(FH);

	my $cmd="
		cat $input_fasta |
		    $REVERSE_COMPLEMENT_FASTA_BIN
			-l $reorient_list >
			$output_fasta
	";

	exec_cmd($cmd);

	# rm reorient_list

}


###############################################################################

sub filter{

	# Calls outside script to filter based on a list.

	my $input_fasta=shift;
	my $keep_hash_ref=shift;
	my $output_fasta=shift;

	open(FH, ">$keep_list") || die "Could not open $keep_list\n";
	foreach my $id(sort keys %{$keep_hash_ref}){
		print FH "$id\n";
	}
	close(FH);

	my $cmd="
		cat $input_fasta |
			$EXTRACT_RECORD_FASTA_BIN
			-l $keep_list >
			$output_fasta
	";

	exec_cmd($cmd);

	# rm keep_list 

}

###############################################################################

sub same_ori{
	# Return 1 if both orientations are the same

	my $q_start=shift;
	my $q_end=shift;
	my $s_start=shift;
	my $s_end=shift;

	# Check query orientation
	my $q_ori;
	if($q_start<$q_end){
		$q_ori=1;
	}else{
		$q_ori=-1;
	}

	# Check subject orientation
	my $s_ori;
	if($s_start<$s_end){
		$s_ori=1;
	}else{
		$s_ori=-1;
	}

	# See if they are in the same orientation
	if($q_ori == $s_ori){
		return(1);
	}else{
		return(0);
	}
}

###############################################################################
