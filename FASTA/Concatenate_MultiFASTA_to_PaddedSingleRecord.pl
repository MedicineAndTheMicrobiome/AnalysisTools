#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_p $opt_o);

my $PADDING=60;

getopts("p:o:");
my $usage = "usage: 
$0 
	-o <ouput filename root>	
	[-p <padding, default=$PADDING>]
	
	This script will read in a multifasta file from STDIN, then generate
	a single fasta by concatenating the records from the multifasta
	file together, with the amount of padding specified. A coordinates
	file will be generated which will specify the ranges of where
	the actual sequences are in the concatenated sequence.  The order
	by which the sequence will be concatenated will be the same as
	the order in which the sequence records were read in.

	The output filename root will also be used as the defline sequence
	id in the output fasta file.

	This script will generate a fasta file and a list of regions that
	real sequence actually exists in.  The file will have the extension:
	.seq_coords.tsv

	Examples:

		$0 -o MyOutput < seq1.fasta

		cat seq1.fasta seq2.fasta | $0 -o MyOutput


	
";

if(!(
	defined($opt_o)
)){
	die $usage;
}

my $output_fname_root = $opt_o;
my $padding = $PADDING;

if(defined($opt_p)){
	$padding = $opt_p;
}
 
###############################################################################

print STDERR "Processing FASTA file...\n\n";

my $total_sequence="";
my @sequence_ranges;
my $cur_begin=0;
my $cur_end=0;

my ($defline, $prev_defline, $sequence);

while(<STDIN>){
	chomp;
	
	if(/^>/){
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

print STDERR "\n";
#print STDERR "$total_sequence\n";

output_fasta_record($output_fname_root, $total_sequence, 80);
output_sequence_coords($output_fname_root, \@sequence_ranges);

print STDERR "\nCompleted.\n";

###############################################################################

sub process_record{
	my $defline=shift;
	my $sequence=shift;
	
	my $seq_len=length($sequence);

	print STDERR "\tAdding: $seq_len residues.\n";	

	if($total_sequence eq ""){
		$cur_begin=0;
		$cur_end=$seq_len;
	}else{
		$sequence= ("N" x $padding) . $sequence;
		$cur_begin=$cur_end+$padding;
		$cur_end=$cur_begin+$seq_len;
	}

	my $seq_id=$defline;
	$seq_id=~s/^>//;
	($seq_id)=split /\s+/, $seq_id;
	push @sequence_ranges, "$seq_id\t$cur_begin\t$cur_end";
	$total_sequence.=$sequence;
}

sub output_fasta_record{
	my $defline = shift;
	my $sequence = shift;
	my $width = shift;

	open(FH, ">$defline\.fasta") || die "Could not open $defline\.fasta\n";

	print FH ">$defline\n";
	my $length=length($sequence);
	print STDERR "Output Sequence Length: $length\n";
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print FH substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

sub output_sequence_coords{
	my $fname=shift;
	my $seq_ranges_ref=shift;

	open(FH, ">$fname\.seq_coords.tsv") || die "Could not open $fname\.seq_coords.tsv\n";

	foreach my $rec (@{$seq_ranges_ref}){
		print FH "$rec\n";		
	}

	close(FH);

}

#------------------------------------------------------------------------------
