#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_m $opt_i $opt_o);

getopts("m:i:o:");
my $usage = "usage: 
$0 
	-i <FASTA input file>
	-o <FASTA output file, with new names>
	-m <Mapping file, old_id -> new_id>

	Reads in a fasta file and then each id gets renamed according to the
	mapping file.  If a new id is not found, then the old id is kept.

";

if(
	!defined($opt_i) ||
	!defined($opt_o) ||
	!defined($opt_m)
){
	die $usage;
}

my $input_fasta=$opt_i;
my $output_fasta=$opt_o;
my $mapping_file=$opt_m;

###############################################################################

my %mapping_hash;

open(MP, "<$mapping_file") || die "Could not open Mapping file: $mapping_file\n";

my $count=0;
while(<MP>){
	chomp;
	my ($old_id, $new_id)=split /\t/, $_;
	$mapping_hash{$old_id}=$new_id;
	$count++;
}

close(MP);

print STDERR "$count mappings loaded.\n";

###############################################################################

open(IN, "<$input_fasta") || die "Could not open input FASTA: $input_fasta\n";
open(OUT, ">$output_fasta") || die "Could not open output FASTA: $output_fasta\n";

print STDERR "Processing FASTA file...\n";

my $mapped=0;
my ($defline, $prev_defline, $sequence);
while(<IN>){
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

print STDERR "$mapped records renamed.\n";
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;

	$defline=~s/^>//;
	my @defarr=split / /, $defline;

	if(defined($mapping_hash{$defarr[0]})){
		$defarr[0]=$mapping_hash{$defarr[0]};
		$mapped++;
	}

	$defline=join / /, @defarr;

	print OUT ">$defline\n";

	my $length=length($sequence);
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print OUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------
