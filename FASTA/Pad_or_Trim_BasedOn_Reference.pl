#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use Sys::Hostname;
use File::Basename;
use FindBin ();
use vars qw($opt_i $opt_r $opt_o $opt_f $opt_t $opt_g);

my $CLUSTALW2_BIN="/usr/local/projects/NTD/SVN/DAS/ViralDataAnalysis/DivergenceMonitor/altnerates/clustalw-2.0.11-linux-i386-libcppstatic/clustalw2";
my $BLASTALL_BIN=`which blastall`;

getopts("i:r:o:f:t:g");
my $usage = "usage: 
$0 
	-i <input fasta>
	-r <reference fasta, this should already be formatdb'd>
	-o <output fasta>

	[-f <number of bases for five prime reference override, eg 12>]
	[-t <number of bases for three prime reference override, eg 13>]

	[-g <ignore errors>]

	This program will go through all the input sequences, compare them against the 
	reference, in order to determine which sequence is most like it.  Then it will
	use clustalw to perform a global alignment and then trim off extra sequence or
	append on extra sequence, depending on what is needed.

	If the -f or -t option is used, then those are the number of bases that wil be
	over ridden by the reference sequence.  These options will fail if there are 
	gaps (indels) in those 5' or 3' ends of the sequence.
	
	This program depends on:

		$BLASTALL_BIN
	and
		$CLUSTALW2_BIN

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
my $five_prime_override=defined($opt_f)?$opt_f:0;
my $three_prime_override=defined($opt_t)?$opt_t:0;
my $ignore_errors=defined($opt_g);

# Echo parameters
print STDERR "Input Fasta: $input_fasta\n";
print STDERR "Output Fasta: $output_fasta\n";
print STDERR "Reference: $reference_fasta\n";

if(defined($five_prime_override)){
	print STDERR "5' override length: $five_prime_override\n";
}
if(defined($three_prime_override)){
	print STDERR "3' override length: $three_prime_override\n";
}

# Get unique temporary file name
my $hostname=hostname();
my $pid=$$;
my ($progname)=fileparse($0);
my $tmp_root="$output_fasta\.$progname\.$hostname\.$pid";

# Set names of intermediate results
my $blast_results_file="$tmp_root\.blastout";
my $aln_results_file="$tmp_root\.alnout";
my $align_input_file="$tmp_root\.alnin";
my $dnd_output_file="$tmp_root\.dnd";

###############################################################################

# Load the target and reference sequences, so we can later break them apart for clustalw
my $input_sequences_hash_ref=load_sequences($input_fasta);
my $reference_sequences_hash_ref=load_sequences($reference_fasta);

# Run blast and try to find best reference sequence for each target sequence
run_blast($input_fasta, $reference_fasta, $blast_results_file);
my ($best_hit_subj_hash, $best_hit_ori_hash)=pick_best_hit($blast_results_file);

open(OUT_FH, ">$output_fasta") || die "Could not open $output_fasta\n";

foreach my $seq_id (sort keys %{$input_sequences_hash_ref}){

	# Build multifasta for clustalw input
	open(FH, ">$align_input_file")  ||  die "Could not open $align_input_file\n";

	# Insert target sequence
	print FH ">$seq_id\n";
	print FH "${$input_sequences_hash_ref}{$seq_id}\n";

	# Make sure a reference sequence is found for the target
	if(!defined(${$best_hit_subj_hash}{$seq_id})){
		die "Could not find hit for $seq_id\n";
	}

	# Insert reference sequence, in proper orientation
	my $best_ref_id=(${$best_hit_subj_hash}{$seq_id});
	my $ref_str=${$reference_sequences_hash_ref}{$best_ref_id};
	if(${$best_hit_ori_hash}{$best_ref_id}==-1){
		$ref_str=revcomp($ref_str);
	}
	print FH ">$best_ref_id\n";
	print FH "$ref_str\n";

	close(FH);

	# Execute alignment
	exec_cmd("$CLUSTALW2_BIN -infile=$align_input_file -outfile=$aln_results_file");
	
	# Parse alignment results
	my $seq_aln_hash_ref=parse_aln_file($seq_id, $best_ref_id, $aln_results_file);

	# Trim or pad gaps
	my $fixed_seq=fixed_gaps(${$seq_aln_hash_ref}{$seq_id}, ${$seq_aln_hash_ref}{$best_ref_id});

	# Override ends with reference
	$fixed_seq=override_sequence($fixed_seq, $five_prime_override, $three_prime_override, ${$seq_aln_hash_ref}{$best_ref_id});
		
	# Remove internal gaps from fixed sequence
	$fixed_seq=~s/-*//g;

	# Write out fixed sequences
	print OUT_FH ">$seq_id\n";
	print OUT_FH "$fixed_seq\n";

}

close(OUT_FH);

# Clean up temporary files
exec_cmd("rm $blast_results_file $aln_results_file $align_input_file $dnd_output_file");

print STDERR "Done.\n\n";


###############################################################################

sub exec_cmd{
	my $cmd=shift;
	$cmd=~s/[\t\n]+/ /g;
	print STDERR "Going to execute: $cmd\n";
	system($cmd);
}

###############################################################################

sub override_sequence{
	my $target_seq=shift;
	my $fiveprime_bases_to_override=shift;
	my $threeprime_bases_to_override=shift;
	my $reference_align_sequence=shift;
	my $overridden_sequence;

	#print "Target: $target_seq\n";
	#print "5' bases: $fiveprime_bases_to_override\n";
	#print "3' bases: $threeprime_bases_to_override\n";

	# remove gaps from sequence ends b/c target should be trim/padded already
	$reference_align_sequence=~s/^-*//;
	$reference_align_sequence=~s/-*$//;
	#print "Reference: $reference_align_sequence\n";

	my $target_seq_length=length($target_seq);
	my $reference_seq_length=length($reference_align_sequence);

	# Get sequences to replace on target
	my $fivep_target_seq=substr($target_seq, 0, $fiveprime_bases_to_override);
	my $threep_target_seq=substr($target_seq, $target_seq_length-$threeprime_bases_to_override, $threeprime_bases_to_override);

	# Get sequences to replace with on reference
	my $fivep_reference_seq=substr($reference_align_sequence, 0, $fiveprime_bases_to_override);
	my $threep_reference_seq=substr($reference_align_sequence, 
		$reference_seq_length-$threeprime_bases_to_override, $threeprime_bases_to_override);
	
	# throw error if there are indels within the override regions, else go ahead and override

	# Work on 5' End
	if($fivep_target_seq=~/-/ || $fivep_reference_seq=~/-/){
		print STDERR "\n";
		print STDERR "****************************************************************\n";
		print STDERR "*     5' target has gaps.  Reference override failed.          *\n";
		print STDERR "****************************************************************\n";
		print STDERR "\n";
		if(!$ignore_errors){
			die "To ignore errors, set ignore errors option flag.\n";
		}
	}else{
		substr($target_seq, 0, $fiveprime_bases_to_override)=$fivep_reference_seq;
	}
	
	# Work on 3' End
	if($threep_target_seq=~/-/ || $threep_reference_seq=~/-/){
		print STDERR "\n";
		print STDERR "****************************************************************\n";
		print STDERR "*     3' target has gaps.  Reference override failed.          *\n";
		print STDERR "****************************************************************\n";
		print STDERR "\n";
		if(!$ignore_errors){
			die "To ignore errors, set ignore errors option flag.\n";
		}
	}else{
		substr($target_seq, $target_seq_length-$threeprime_bases_to_override, $threeprime_bases_to_override)=$threep_reference_seq;
	}

	return($target_seq);
}

###############################################################################

sub load_sequences{
	my $file=shift;
	my %sequence_hash;	

	print STDERR "Loading $file\n";
	open(FH, "<$file") || die "Could not open $file for reading.\n";

	my $sequence="";
	my $defline;
	my $prev_defline;
	while(<FH>){
		chomp;

		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, \$sequence, \%sequence_hash);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	process_record($prev_defline, \$sequence, \%sequence_hash);

	close(FH);

	return(\%sequence_hash);
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

	${$length_hash_ref}{$id}=${$sequence_ref};
	print STDERR ".";
}

###############################################################################

sub run_blast{
	my $input_file=shift;
	my $reference_fasta=shift;
	my $results_file=shift;

	my $cmd="
		$BLASTALL_BIN 
			-p blastn
			-i $input_file 
			-d $reference_fasta
			-m 8 
			-e 1e-100
			-o $results_file
	";
	exec_cmd($cmd);

}

###############################################################################

sub pick_best_hit{
	my $results_file=shift;

	open(FH, "<$results_file") || die "Could not open $results_file for reading.\n";

	my %best_hit_hash;
	my %best_hit_subject_hash;
	my %best_hit_ori_hash;

	while(<FH>){
		chomp;
		
		# Parse the blast results
		my 
			($query_id,$subject_id,$perc_identity,$alignment_length,
			$mismatches,$gap_openings,$q_start,$q_end,$s_start,$s_end,
			$evalue,$bit_score)=split /\t/, $_;
				
		if(!defined($best_hit_hash{$query_id})){
			# If sequence not seen yet, take first as best/only hit
			$best_hit_hash{$query_id}=$bit_score;
			$best_hit_subject_hash{$query_id}=$subject_id;
			$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
		}else{
			# Replace best if we find something with a higher bit score
			my $prior_score=$best_hit_hash{$query_id};
			if($prior_score<$bit_score){
				$best_hit_hash{$query_id}=$bit_score;
				$best_hit_subject_hash{$query_id}=$subject_id;
				$best_hit_ori_hash{$query_id}=same_ori($q_start,$q_end,$s_start,$s_end);
			}
		}
	}

	close(FH);

	# Return a hash of query id's that are in different orientation with reference/subject
	return(\%best_hit_subject_hash, \%best_hit_ori_hash);
	
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

sub revcomp{
    my $sequence = shift;
    $sequence =~ tr/atugcyrswkmbdhvnxATUGCYRSWKMBDHVNX/taacgryswmkvhdbnxTAACGRYSWMKVHDBNX/;
    my $revcomp = "";
    for (my $j=length($sequence)-1;$j>-1;$j--) {
        $revcomp .= substr($sequence,$j,1);
    }
    return $revcomp;
}

###############################################################################

sub parse_aln_file{
	my $target_seq_id=shift;
	my $reference_seq_id=shift;
	my $aln_results_file=shift;
	my $fixed_seq;

	my %algn_str_hash;

	open(FH, "<$aln_results_file") || die "Could not open $aln_results_file\n";
	while(<FH>){
		chomp;
		my ($seq_id, $algn_str)=split /\s+/, $_;
		$algn_str_hash{$seq_id}.=$algn_str;
	}
	close(FH);

	return(\%algn_str_hash);

}

###############################################################################

sub fixed_gaps{

	my $tar_seq=shift;
	my $ref_seq=shift;

	#print STDERR "\n";	
	#print STDERR "$tar_seq\n\n";
	#print STDERR "$ref_seq\n";

	my $tar_begin_gaps=num_begin_gaps($tar_seq);
	my $ref_begin_gaps=num_begin_gaps($ref_seq);

	print STDERR "Beginning: target gaps = $tar_begin_gaps, reference gaps = $ref_begin_gaps\n";

	my $tar_end_gaps=num_end_gaps($tar_seq);
	my $ref_end_gaps=num_end_gaps($ref_seq);

	print STDERR "Ending: target gaps = $tar_end_gaps, reference gaps = $ref_end_gaps\n";

	$tar_seq=~s/^-*//;
	$tar_seq=~s/-*$//;
	$ref_seq=~s/^-*//;
	$ref_seq=~s/-*$//;

	my $clean_tar_len=length($tar_seq);
	my $clean_ref_len=length($ref_seq);

	my $fixed_seq;

	# Fix beginning gaps
	if($tar_begin_gaps>0){
		my $front_pad=substr($ref_seq, 0, $tar_begin_gaps);
		$fixed_seq=$front_pad . $tar_seq;
	}elsif($ref_begin_gaps>0){
		$fixed_seq=$tar_seq;
		substr($fixed_seq, 0, $ref_begin_gaps)="";
	}else{
		$fixed_seq=$tar_seq;
	}

	# Fix ending gaps
	if($tar_end_gaps>0){
		$fixed_seq.=substr($ref_seq, $clean_ref_len-$tar_end_gaps, $tar_end_gaps);
	}elsif($ref_end_gaps>0){
		my $fixed_seq_len=length($fixed_seq);
		substr($fixed_seq, $fixed_seq_len-$ref_end_gaps, $ref_end_gaps)="";
	}

	return($fixed_seq);
}

#------------------------------------------------------------------------------

sub num_end_gaps{
	my $seq=shift;
	if($seq=~/(-*)$/){
		return(length($1));
	}else{
		return(0);
	}
}

#------------------------------------------------------------------------------

sub num_begin_gaps{
	my $seq=shift;
	if($seq=~/^(-*)/){
		return(length($1));
	}else{
		return(0);
	}
}

#------------------------------------------------------------------------------

