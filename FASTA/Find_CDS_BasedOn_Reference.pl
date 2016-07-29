#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use Sys::Hostname;
use File::Basename;
use FindBin ();
use vars qw($opt_i $opt_f $opt_t $opt_o);

my $BLASTALL_BIN=`which blastall`;

getopts("i:f:t:o:");
my $usage = "usage: 
$0 
	-i <input fasta>
	-f <five prime reference fasta, this should already be formatdb'd>
	-t <three prime reference fasta, this should already be formatdb'd>
	-o <output CDS positions>

	This program will compute the begin and end of where the reference
	AA sequences are found in the input nucleotide fasta.

	The best hit from EACH of the 2 fasta files will be used to compute the outer
	bounds of the CDS.  It's not elegant and requires curation of the
	two reference fastas, but it should be fairly robust.

	Output looks like:

	    <seq_id>\\t<dist_from_begin>\\t<dist_from_end>\\t<seq_length\\n

	This program depends on:

		$BLASTALL_BIN
";

if(!(
	defined($opt_i) ||
	defined($opt_f) ||
	defined($opt_t) ||
	defined($opt_o))){
	die $usage;
}

my $input_fasta=$opt_i;
my $reference_fasta_5p=$opt_f;
my $reference_fasta_3p=$opt_t;
my $output_cds=$opt_o;

# Echo parameters
print STDERR "\n";
print STDERR "Input Fasta: $input_fasta\n";
print STDERR "Output Fasta: $output_cds\n";
print STDERR "Reference 5prime: $reference_fasta_5p\n";
print STDERR "Reference 3prime: $reference_fasta_3p\n";
print STDERR "\n";

# Get unique temporary file name
my $hostname=hostname();
my $pid=$$;
my ($progname)=fileparse($0);
my $tmp_root="$output_cds\.$progname\.$hostname\.$pid";

# Set names of intermediate results
my $blast_results_file="$tmp_root\.blastout";

###############################################################################

# Get sequence lengths
my %input_len_hash;
my %ref_len_hash;
get_sequence_lengths(\%input_len_hash, $input_fasta);
get_sequence_lengths(\%ref_len_hash, $reference_fasta_5p);
get_sequence_lengths(\%ref_len_hash, $reference_fasta_3p);

my @reference_files=($reference_fasta_5p, $reference_fasta_3p);
my %results_hash;

foreach my $ref_file(@reference_files){

	# Run blast
	run_blast($input_fasta, $ref_file, $blast_results_file);

	# Get best hit
	my ($best_hit_ori_hash_ref, $best_hit_sbj_aln_hash_ref,
		 $best_hit_ref_aln_hash_ref, $best_hit_ref_hash_ref)=parse_blast_for_best_hit($blast_results_file);

	# Get coordinates for CDS
	get_coordinates(
		\%results_hash,
		\%input_len_hash, \%ref_len_hash, 
		$best_hit_ref_hash_ref,
		$best_hit_ori_hash_ref, 
		$best_hit_ref_aln_hash_ref, $best_hit_sbj_aln_hash_ref);

	print STDERR "\n";

}

# Take outer most of coordinates
my $outer_coords_hash_ref=get_outer_coords(\%results_hash);

# Output results 
open(OUTFH, ">$output_cds") || die "Could not open $output_cds for writing.\n";

# Compute distance from ends and output
foreach my $query_id(sort keys %{$outer_coords_hash_ref}){

	my ($begin, $end)=split /#/, ${$outer_coords_hash_ref}{$query_id};
	my $query_length=$input_len_hash{$query_id};

	my $front_dist=$begin;
	my $end_dist=$query_length-$end;

	print OUTFH "$query_id\t$front_dist\t$end_dist\t$query_length\n";
}

close(OUTFH);

print STDERR "Done.\n";

###############################################################################

sub exec_cmd{
	my $cmd=shift;
	$cmd=~s/[\t\n]+/ /g;
	print STDERR "Going to execute: $cmd\n";
	system($cmd);
}

###############################################################################

sub get_outer_coords{
	my $results_hash=shift;
	my %outer_coords;
			
	foreach my $query_id(sort keys %{$results_hash}){

		my @coords_arr=@{${$results_hash}{$query_id}};

		my ($min, $max)=split /#/, $coords_arr[0];

		foreach my $coord(@coords_arr){
			my ($begin, $end)=split /#/, $coord;
			if($begin<$min){
				$min=$begin;
			}
			if($end>$max){
				$max=$end;
			}	
		}

		$outer_coords{$query_id}="$min#$max";
	}

	return(\%outer_coords);
}

###############################################################################

sub get_coordinates{
	my $results_hash_ref=shift;
	my $qry_len_hash_ref=shift;
	my $ref_len_hash_ref=shift;
	my $best_hit_ref_hash_ref=shift;
	my $best_hit_ori_hash_ref=shift;
	my $best_hit_ref_aln_hash_ref=shift;
	my $best_hit_qry_aln_hash_ref=shift;

	foreach my $query_id(sort keys %{$qry_len_hash_ref}){

		print STDERR "\tWorking on $query_id\n";

		my $reference_id=${$best_hit_ref_hash_ref}{$query_id};
		my $ori=${$best_hit_ori_hash_ref}{$query_id};
		my ($rbegin, $rend)=split /#/, ${$best_hit_ref_aln_hash_ref}{$query_id};
		my ($qbegin, $qend)=split /#/, ${$best_hit_qry_aln_hash_ref}{$query_id};
		my $ref_length=${$ref_len_hash_ref}{$reference_id};
		my $qry_length=${$qry_len_hash_ref}{$query_id};


		if(!defined($reference_id)){
			print STDERR "\tCould not find a reference hit for $query_id\n";
			next;
		}

		if($ori != 1){
			die "Orientation between best hit and sequence are not consistent for $query_id\n";
			next;
		}

		if($qbegin > $qend){
			die "Orientation of query sequence is not in forward direction for $query_id\n";
			next;
		}

		#print STDERR "$reference_id\n";
		#print STDERR "$ori\n";
		#print STDERR "r: $rbegin - $rend\n";
		#print STDERR "q: $qbegin - $qend\n";
		#print STDERR "q: $qry_length ref: $ref_length\n";
	
		# Pad to same length as reference
		my $upstream_pad=$rbegin*3;
		my $downstream_pad=($ref_length-$rend)*3;
		
		# Add stop codon
		$downstream_pad+=3;

		$qbegin-=$upstream_pad;
		$qend+=$downstream_pad;

		#print STDERR "$qbegin / $qend\n";

		#my $cds_upend=$qbegin;
		#my $cds_downend=$qry_length - $qend;
		#print STDERR "$cds_upend :: $cds_downend\n";

		push @{${$results_hash_ref}{$query_id}}, "$qbegin#$qend";

	}
	
}

###############################################################################

sub get_sequence_lengths{
	my $length_hash_ref=shift;	
	my $file=shift;

	open(FH, "<$file") || die "Could not open $file for reading.\n";

	my $sequence="";
	my $defline;
	my $prev_defline;
	while(<FH>){
		chomp;

		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, \$sequence, $length_hash_ref);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	process_record($prev_defline, \$sequence, $length_hash_ref);

	close(FH);
}

#------------------------------------------------------------------------------

sub process_record{
	my $defline=shift;
	my $sequence_ref=shift;
	my $length_hash_ref=shift;
	
	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Error parsing deflined: $defline\n";
	}

	${$length_hash_ref}{$id}=length(${$sequence_ref});
}

###############################################################################

sub run_blast{
	my $input_file=shift;
	my $reference_fasta=shift;
	my $results_file=shift;

	my $cmd="
		$BLASTALL_BIN 
			-p blastx
			-i $input_file 
			-d $reference_fasta
			-m 8 
			-e 1e-10
			-o $results_file
	";
	exec_cmd($cmd);

}

###############################################################################

sub parse_blast_for_best_hit{
	my $results_file=shift;
	my $ref_len_hash_ref=shift;
	my $filter_cutoff=shift;

	open(FH, "<$results_file") || die "Could not open $results_file for reading.\n";

	my %best_hit_hash;
	my %best_hit_ori_hash;
	my %best_hit_ref_aln_hash;
	my %best_hit_sbj_aln_hash;
	my %best_hit_ref_hash;

	while(<FH>){
		chomp;
		
		# Parse the blast results
		my 
			($query_id,$subject_id,$perc_identity,$alignment_length,
			$mismatches,$gap_openings,$q_start,$q_end,$s_start,$s_end,
			$evalue,$bit_score)=split /\t/, $_;

			# Convert to 0-space based coordinates
			if($q_start <= $q_end){
				$q_start--;
			}else{
				$q_end--;
			}
				
			if($s_start <= $s_end){
				$s_start--;
			}else{
				$s_end--;
			}
				
		if(!defined($best_hit_hash{$query_id})){
			# If sequence not seen yet, take first as best/only hit
			$best_hit_hash{$query_id}=$bit_score;
			$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
			$best_hit_ref_hash{$query_id}=$subject_id;
			$best_hit_ref_aln_hash{$query_id}="$s_start#$s_end";
			$best_hit_sbj_aln_hash{$query_id}="$q_start#$q_end";
		}else{
			# Replace best if we find something with a higher bit score
			my $prior_score=$best_hit_hash{$query_id};
			if($prior_score<$bit_score){
				$best_hit_hash{$query_id}=$bit_score;
				$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
				$best_hit_ref_hash{$query_id}=$subject_id;
				$best_hit_ref_aln_hash{$query_id}="$s_start#$s_end";
				$best_hit_sbj_aln_hash{$query_id}="$q_start#$q_end";
			}
		
		}

	}

	close(FH);

	# Return a hash of query id's that are in different orientation with reference/subject
	return(\%best_hit_ori_hash, \%best_hit_sbj_aln_hash, \%best_hit_ref_aln_hash, \%best_hit_ref_hash);
	
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
