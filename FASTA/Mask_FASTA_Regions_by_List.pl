#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_l $opt_M $opt_w $opt_b $opt_e $opt_r);

my $WIDTH=60;
my $DEFAULT_BEGIN=1;
my $DEFAULT_END=2;

getopts("f:l:M:w:b:e:r");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-l <list of masking coordinates>
	[-M <mask char/quality value>]
	[-w <width, default $WIDTH>]
	[-r Assume input is in 1-residue based coordinates]

	Use alternate trim columns, counting from 0.
	[-b <clear begin column, default=$DEFAULT_BEGIN>]
	[-e <clear end column, default=$DEFAULT_END>]

	Given a list of regions to mask, ie.

	<seq_id>\\t<begin>\\t<end>\\n

	the script will apply all the masking and write 
	the sequence to STDOUT.

";

if(!(
	defined($opt_f) && 
	defined($opt_l))){
	die $usage;
}

my $mask_char="N";
if(defined($opt_M)){
	$mask_char=$opt_M;
}

my $width=$WIDTH;
if(defined($opt_w)){
	$width=$opt_w;
}

my $begin_col=$DEFAULT_BEGIN;
if(defined($opt_b)){
	$begin_col=$opt_b;
}

my $end_col=$DEFAULT_END;
if(defined($opt_e)){
	$end_col=$opt_e;
}

my $residue_based=0;
if(defined($opt_r)){
	$residue_based=1;
	print STDERR "Being told coordinates are in 1-Residue-based coordinates.\n";
}


###############################################################################
# Read in trim list

print STDERR "Loading trim list: $opt_l\n";
open(LIST_FH, "<$opt_l") || die "Could not open $opt_l\n";

my %list_hash;
my $list_length=0;
my $num_regions=0;
while(<LIST_FH>){
	chomp;

	my @in=split /\t/, $_;
	my $id=$in[0];
	my $begin=$in[$begin_col];
	my $end=$in[$end_col];

	push @{$list_hash{$id}}, "$begin#$end";
	$num_regions++;
}

close(LIST_FH);

print STDERR "Done.\n";
print STDERR "$num_regions regions loaded.\n";
my $num_sequences_to_mask=keys %list_hash;
print STDERR "$num_sequences_to_mask targeted for masking.\n";

###############################################################################
# Read in features

my $num_found=0;
my $num_regions_masked=0;
my $num_sequences_masked=0;

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

print STDERR "Processing FASTA file...\n";
my $num_sequences_in_fasta=0;
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
	chomp;
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence, $residue_based);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

close(FASTA_FH);

print STDERR "$num_sequences_in_fasta sequences in input FASTA\n";
print STDERR "$num_regions_masked of $num_regions targeted regions masked\n";
print STDERR "$num_sequences_masked of $num_sequences_to_mask targeted sequences masked \n";
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse id from defline: $defline\n";
	}

	# Remove leading and trailing spaces
	$sequence=~s/^\s+//;
	$sequence=~s/\s+$//;

	if(defined($list_hash{$id})){
		foreach my $mask_region(@{$list_hash{$id}}){
			my ($begin, $end)=split /#/, $mask_region;
			$sequence=mask(\$sequence, $begin, $end, $residue_based);
			$num_regions_masked++;
		}
		$num_sequences_masked++;
	}
	
	print STDOUT "$defline\n";
	dump_sequence(\$sequence);

	$num_sequences_in_fasta++;
}

#------------------------------------------------------------------------------

sub mask{
	my $seq_ref=shift;
	my $begin=shift;
	my $end=shift;
	my $residue_based=shift;

	my $seq=${$seq_ref};

	if($begin>$end){
		($end,$begin)=($begin,$end);
	}
	if($residue_based){
		$begin--;
	}

	#print STDERR "Masking from $begin - $end\n";

	my $mask_length=$end-$begin;
	my $substitution_string=($mask_char x $mask_length);
	substr($seq, $begin, $mask_length)=$substitution_string;
	return($seq);

}
	
sub dump_sequence{
	my $seq_ref=shift;

	my $length=length(${$seq_ref});
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr(${$seq_ref}, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
		
}
