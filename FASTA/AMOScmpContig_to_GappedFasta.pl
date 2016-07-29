#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_i $opt_o);

getopts("i:o:");
my $usage = "usage: 
$0 
	-i <AMOScomp Contig File>
	-o <Out name root>

	This program will read in an AMOScomp contig file and convert it
	into a gapped FASTA file.  There should be only one contig (with
	many reads) in the AMOScomp contig file.  Since AMOScomp performs
	trimming, a separate file will be generated which contains the
	range of the read that was actually aligned.

	<read id>\\t<used_begin>\\t<used_end>\\t<orientation>\\n

	Output will be saved to:

	<Out name root>\.gap_fasta
	<Out name root>\.reads_info

";

if(
	!(defined($opt_i)) ||
	!(defined($opt_o))
){
	die $usage;
}

my $amoscomp_in=$opt_i;
my $gappedfasta_out="$opt_o\.gap_fasta";
my $readsinfo_out="$opt_o\.reads_info";

my $gapped_consensus_length;
my @id_arr;
my %align_hash;
my %sequence_hash;
my %read_used_range_hash;

###############################################################################

print STDERR "Processing AMOScomp file...\n";

open(AMOS_IN, "<$amoscomp_in") || die "Could not open $amoscomp_in\n";
open(GAPPED_FASTA_OUT, ">$gappedfasta_out") || die "Could not open $gappedfasta_out\n";
open(READS_INFO_OUT, ">$readsinfo_out") || die "Could not open $readsinfo_out\n";

my ($defline, $prev_defline, $sequence);
while(<AMOS_IN>){
	chomp;
	
	if(/^#/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

print STDERR "Gapped consensus length: $gapped_consensus_length\n";

# Check to make sure sequences are the same length


my $max_length=0;
my %padded_seq_hash;
for(my $i=0; $i<=$#id_arr; $i++){
	my $id=$id_arr[$i];
	my ($begin, $end)=split /#/, $align_hash{$id};

	# Determine how much to pad 5' and 3' end
	my $padded_seq = (('-' x $begin) . $sequence_hash{$id} . ('-' x ($gapped_consensus_length-$end)));
	my $padded_seq_len=length($padded_seq);
	$padded_seq_hash{$id}=$padded_seq;

	# Keep track of longest sequence
	if($max_length<$padded_seq_len){
		$max_length=$padded_seq_len;
	}
}

for(my $i=0; $i<=$#id_arr; $i++){
	my $id=$id_arr[$i];

	# Determine if we need to do extra padding
	my $padded_seq=$padded_seq_hash{$id};
	my $diff_from_max_length=$max_length-length($padded_seq);
	if($diff_from_max_length>0){
		$padded_seq.=("-" x $diff_from_max_length);
		print STDERR "WARNING: Padding 3' end with -.\n";
	}

	# Output padded sequences
	output_fasta_record(\*GAPPED_FASTA_OUT, ">$id", $padded_seq);
	print READS_INFO_OUT "$id\t$read_used_range_hash{$id}\n";
}

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	my $gapped_read_align_begin;
	my $gapped_read_align_length;
	my $gapped_read_align_end;
	my $read_begin;
	my $read_end;
	my $read_ori;

	if($defline=~/^##/){
		# Skip consensus line, but
		# get the gapped consensus length

		if($defline=~/ (\d+) bases/){
			$gapped_consensus_length=$1;
		}else{
			die "Could not parse consensus record for gapped consensus length: '$defline'\n";
		}
		
		return;
	}

	if($defline=~/^#([_A-Z0-9]+)/){
		$id=$1;	
		#print STDERR "$id\n";
		if($defline=~/\((\d+)\)/){
			$gapped_read_align_begin=$1;
		}else{
			print STDERR "Could not parse alignment extents: '$defline'\n";
		}

		if($defline=~/ (\d+) bases/){
			$gapped_read_align_length=$1;
			$gapped_read_align_end=$gapped_read_align_begin+$gapped_read_align_length;
		}else{
			die "Could not parse consensus record for gapped consensus length: $defline\n";
		}

		if($defline=~/ {(\d+) (\d+)} /){
			$read_begin=$1;
			$read_end=$2;
			$read_ori;
			
			# Convert to begin/end/orientation
			if($read_begin>$read_end){
				($read_begin, $read_end)=($read_end, $read_begin);
				$read_ori=-1
			}else{
				$read_ori=1;
			}		
			$read_begin--;
		}		

	}else{
		print STDERR "Could not parse '$defline'\n";
	}

	# Store sequence and alignment information
	push @id_arr, $id;
	$align_hash{$id}="$gapped_read_align_begin#$gapped_read_align_end";
	$sequence_hash{$id}=$sequence;
	$read_used_range_hash{$id}="$read_begin\t$read_end\t$read_ori";

	return;
}

###############################################################################

sub output_fasta_record{
	my $fh=shift;
	my $defline=shift;
	my $sequence=shift;

	print {$fh} "$defline\n";

	my $length=length($sequence);
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print {$fh} substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------
