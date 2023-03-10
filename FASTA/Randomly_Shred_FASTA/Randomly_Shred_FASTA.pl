#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_n $opt_l $opt_o);

getopts("f:n:l:o:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-n <Num Shreds/Fragments>
	-l <length>
	-o <output file name>

	This program will take the input fasta file and randomly generate
	a set of reads with the specified length.

";

if(!(
	defined($opt_f) && 
	defined($opt_n) && 
	defined($opt_l) && 
	defined($opt_o))){
	die $usage;
}

###############################################################################

my $fasta=$opt_f;
my $num_shreds=$opt_n;
my $shred_length=$opt_l;
my $output_root=$opt_o;

print STDERR "Input fasta: $fasta\n";
print STDERR "Number of shreds: $num_shreds\n";
print STDERR "Shred Length: $shred_length\n";
print STDERR "Output root: $output_root\n";

###############################################################################

open(OUT_FH, ">$output_root") || die "Could not open $output_root\n";

open(FASTA_FH, "<$fasta") || die "Could not open $fasta\n";

my $offset=0;
print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
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

close(FASTA_FH);

print STDERR "Completed.\n";

###############################################################################
###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $seq_len=length($sequence);

	print STDERR "Processing:\n";
	print STDERR "$defline\n";
	print STDERR "Sequence Length: $seq_len\n";

	my @splits=split /\s+/, $defline;
	my $clean_id=$splits[0];

	for(my $i=0; $i<$num_shreds; $i++){

		# Find fragment to extract as shred
		my $rand_pos=int(rand()*($seq_len-$shred_length));
		print "Shred Pos: $rand_pos\n";
		my $frag=substr($sequence, $rand_pos, $shred_length);

		# Output Shred	
		my $shred_id=sprintf("%03i", $i+1);
		print OUT_FH "$clean_id\.$shred_id\n";
		my $length=length($frag);
		my $width=70;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print OUT_FH substr($frag, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);

	}

}

#------------------------------------------------------------------------------

