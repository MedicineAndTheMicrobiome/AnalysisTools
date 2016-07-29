#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_a $opt_o $opt_t);
use File::Basename;

my $THRESHOLD=100;

getopts("a:o:t:");
my $usage = "usage: 
$0 
	-a <.ALN file>
	[-t <representation threshold, default $THRESHOLD]
	[-o <output file name, default <.ALN file>.fasta >]

	Takes aln file from clustalw and generates a consensus fasta file.

	This program is similar to the consambig from EMBOSS, except that
	instead N's are only output if there is evidence for all {A,T,G,C}.  If
	N's exist in the input, then they are ignored.  Gaps represented by
	hyphens (-) in the input will be treated as N's.  N's will only be generated
	for output if there is evidence for all 4 nucleotides in the input
	sequence, otherwise the best 2 or 3 based encoding ambiguiting code
	will be chosen.

	The sequence id used in the defline will be the file name of the align
	file, excluding the .aln extension, and excluding the path of the aln file.

	The -t <representation threshold>:
		(100 - threshold)% is the percentage representation that must be
		present in the population for that variation to not be considered
		a sequencing error.  So if you set the -t option to 99, then if 
		an allele occurs with less than a rate of 1%, it will not be
		represented by an ambiguity code.
	
";

if(!defined($opt_a)){
	die $usage;
}

my $aln_file=$opt_a;
my $threshold=defined($opt_t)?$opt_t:$THRESHOLD;
my $output_file=$aln_file;

my $min_repr_threshold=(100.0-$threshold)/100.0;
print STDERR "Min representation threshold: $min_repr_threshold\n";

if(defined($opt_o)){
	$output_file=$opt_o;
}else{
	$output_file=~s/aln$/$threshold.cons.fasta/;
}
print STDERR "Output is going to: $output_file\n";

my ($seqname, $path)=fileparse($aln_file);
$seqname=~s/\.aln//;

###############################################################################

print STDERR "Reading ALN file...\n";

open(ALN_FH, "<$aln_file") || die "Could not open $aln_file\n";

my %seq_hash;
my $clustal_file=0;
while(<ALN_FH>){
	chomp;
	if($_=~/CLUSTAL/ || $_=~/MUSCLE/){
		$clustal_file=1;
		next;	
	}elsif($_=~/\*/){
		next;
	}else{
		my ($id, $sequence)=split /\s+/, $_;
		if($id ne ""){
			$seq_hash{$id}.=uc($sequence);
		}
	}
}
close(ALN_FH);

if($clustal_file!=1){
	print STDERR "Are you sure this is clustal file?\n";
}

print STDERR "done.\n";

###############################################################################

print STDERR "Computing nucleotide profiles...\n";
my @nuc_at_pos;
my $length=-1;
foreach my $id(keys %seq_hash){
	#print ">$id<\n";
	#print "$seq_hash{$id}\n";

	my @nucs=split //, $seq_hash{$id};
	
	# Assume length of consensus is the length of the first sequence, but die if not all
	#  lengths are the same length.
	my $curlength=$#nucs+1;
	if($length==-1){
		$length=$curlength;
	}elsif($length != $curlength){
		die "Error, for some reason the length of alignment for $id is not the same as all the other sequences ($length != $curlength)\n";
	}

	# Compute the profile at each position
	#	Each position is represented in an array
	#	Each nucleotide code is summed up in a hash
	for(my $i=0; $i<$length; $i++){
		${$nuc_at_pos[$i]}{$nucs[$i]}++;
	}
}
print STDERR "done.\n";

###############################################################################

print STDERR "Normalizing profiles and eliminating nucs that don't pass threshold...\n";
for(my $i=0; $i<$length; $i++){

	# Sum up the total number of nucleotides at each position
	my $tot=0;
	foreach my $nuc (keys %{$nuc_at_pos[$i]}){
		$tot+=${$nuc_at_pos[$i]}{$nuc};
	}

	# Eliminated the variation at each position if it doesn't exceed the min_repr_threshold
	foreach my $nuc (keys %{$nuc_at_pos[$i]}){
		${$nuc_at_pos[$i]}{$nuc}/=$tot;
		if( $min_repr_threshold > ${$nuc_at_pos[$i]}{$nuc} ){
			delete ${$nuc_at_pos[$i]}{$nuc};
		}
	}
}

###############################################################################

print STDERR "Building consensus sequence...\n";

my $consensus="";
for(my $i=0; $i<$length; $i++){
	$consensus.=get_amb($nuc_at_pos[$i]);
}

###############################################################################

print STDERR "Writing Consensus: $output_file\n";

output_fasta($output_file, $consensus);

print STDERR "Done!\n";

###############################################################################

sub get_amb{
	my $nuc_hash_ref=shift;
	my $code;

	# The order of these if statements are important.

	if(has($nuc_hash_ref, "A")){
		$code="A";
	}
	if(has($nuc_hash_ref, "C")){
		$code="C";
	}
	if(has($nuc_hash_ref, "G")){
		$code="G";
	}
	if(has($nuc_hash_ref, "T")){
		$code="T";
	}
	if(has($nuc_hash_ref, "AC")){
		$code="M";
	}
	if(has($nuc_hash_ref, "AG")){
		$code="R";
	}
	if(has($nuc_hash_ref, "AT")){
		$code="W";
	}
	if(has($nuc_hash_ref, "CG")){
		$code="S";
	}
	if(has($nuc_hash_ref, "CT")){
		$code="Y";
	}
	if(has($nuc_hash_ref, "GT")){
		$code="K";
	}
	if(has($nuc_hash_ref, "ACG")){
		$code="V";
	}
	if(has($nuc_hash_ref, "ACT")){
		$code="H";
	}
	if(has($nuc_hash_ref, "AGT")){
		$code="D";
	}
	if(has($nuc_hash_ref, "CGT")){
		$code="B";
	}
	if(has($nuc_hash_ref, "GATC")){
		$code="N";
	}
	return($code);
}

###############################################################################

sub has{
	# See if all the members in the string, have a representive in the hash

	my $hash_ref=shift;
	my %hash=%{$hash_ref};	# Hash of nucleotide counts
	my $str=shift;		# String containing representatives of ambiguity code

	# Test all the nucs in the hash to see if they exist in the hash
	my @nucs=split //, $str;
	my $num_hits=0;
	foreach my $nuc(@nucs){
		if(defined($hash{$nuc})){
			$num_hits++;
		}	
	}

	# Only return true if all the members in the string are represented in the hash
	if($num_hits==($#nucs+1)){
		return(1);
	}else{
		return(0);
	}
}

###############################################################################

sub output_fasta{
	my $filename=shift;
	my $sequence=shift;

	open(FH, ">$filename") || die "Could not open $filename\n";

	print FH ">$seqname\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print FH substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

	close(FH);
}

#------------------------------------------------------------------------------
