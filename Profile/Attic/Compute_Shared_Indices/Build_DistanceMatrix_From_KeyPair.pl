#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use vars qw($opt_c $opt_o $opt_s);
use FindBin;
use FileHandle;
getopts("c:o:s");

my $usage = "usage:
$0
	
	[-c <column number to find distance value, column 0 and 1 are the ids for which the distance refers to>]
		(default is 2, the 3rd column.)
	[-o <output distance matrix>]
	[-s (flag to convert from similarity to distance)]

	<List of key-pair value files>
		

	Example Input:

	v35.159733294.zo9_Tongue_dorsum,v35.158398106.zg_Stool,0.0258964082283755
	v35.159733294.zo9_Tongue_dorsum,v35.158398106.zo1_Buccal_mucosa,0.50576088111498
	v35.159733294.zo9_Tongue_dorsum,v35.158398106.zo2_Hard_palate,0.458747734143091
	v35.159733294.zo9_Tongue_dorsum,v35.158398106.zo3_Attached_gingivae,0.341735271245295

";

print STDERR $usage;

my $column=$opt_c;
my $outmat=$opt_o;
my $similarity=defined($opt_s);

if($similarity){
	print STDERR "Planning on converting similarities to distances.\n";
}

if(!defined($column)){
	$column=2;
}
if(!defined($outmat)){
	$outmat="output.mat";
}

my $delim=",";

print STDERR "Working on column $column\n";

my %all_entries_hash;
my %all_keys_hash;

my $infile;
while($infile=shift){
	print STDERR "Working on $infile...\n";

	open(IN_FH, "<$infile") || die "Could not open $infile.\n";

	while(<IN_FH>){
		chomp;
	
		if($_=~/^#/){
			print STDERR "Comment/Header line detected.\n";
			next;
		}

		my @columns=split /$delim/,$_;

		my $key1=$columns[0];
		my $key2=$columns[1];
		my $dist=$columns[$column];

		if($similarity){
			$dist=1-$dist;	
		}

		$all_entries_hash{$key1}{$key2}=$dist;
		
		$all_keys_hash{$key1}=1;
		$all_keys_hash{$key2}=1;
	}

}

open(OUTFH, ">$outmat") || die "Could not open $outmat\n";

my @keys=sort keys %all_keys_hash;


# Print column labels
my $col=0;
foreach my $key(@keys){
	print OUTFH ",$key";
	$col++;
}
print OUTFH "\n";

# Print distances
foreach my $key1 (@keys){
	my $col=0;
	print OUTFH $key1;
	foreach my $key2 (@keys){
		print OUTFH ",$all_entries_hash{$key1}{$key2}";		
		$col++;
	}
	print OUTFH "\n";
}

close(OUTFH);
