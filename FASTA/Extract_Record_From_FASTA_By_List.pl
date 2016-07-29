#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_l);

getopts("f:l:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-l <list of defline identifiers>

	Given a list of identifiers and a FASTA file, this program
	will report the FASTA records that match to STDOUT.

	This is an exact match of the defline identifier, so it
	ignores any defline tags.  

	So if the input identifiers file contains EXAMPLE, then
	deflines that look like:

	>EXAMPLE /tag1=xxx /tag2=xxx

	or

	>EXAMPLE /tag3=xxx

	will be reported, but

	>EXAMPLE1

	will not be reported.

";

if(!(
	defined($opt_f) && 
	defined($opt_l))){
	die $usage;
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(LIST_FH, "<$opt_l") || die "Could not open $opt_l\n";

my %list_hash;
my $list_length=0;
while(<LIST_FH>){
	chomp;
	my ($id)=split /\t/, $_;
	$list_hash{$id}=1;
	$list_length++;
}

###############################################################################
# Read in features

my $num_found=0;

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

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

print STDERR "$num_found out of $list_length found\n";
print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	$defline=~/^>(\S+)/;
	my $id=$1;
	
	if(defined($list_hash{$id})){

		print STDOUT "$defline\n";

		my $length=length($sequence);
		my $width=200;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print STDOUT substr($sequence, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);
		
		$num_found++;
	}
}

#------------------------------------------------------------------------------

