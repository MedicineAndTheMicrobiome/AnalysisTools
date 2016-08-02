#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_m);

getopts("i:m:");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-m <Map file, tab-separated>


	This script will read in a reaction file and then generate
	a bounds file with the following format per line:

	<BiGG Metabolite ID> \\t <lowbnd> \\t <uppbnd> \\t <description> \\n

	The purpose of this script is to extract the bounds from the
	reaction file into a bounds file with the 'universal' BiGG ID
	mapping.

";

if(!defined($opt_i) || !defined($opt_m)){
	die $usage;
}

my $InputFile=$opt_i;
my $MapFile=$opt_m;
my $OutputFile=$InputFile;

if($OutputFile=~/\.bnds\.tsv$/){
	$OutputFile=~s/\.bnds\.tsv$/.BiGG.bounds.tsv/;
}elsif($OutputFile=~/\.tsv$/){
	$OutputFile=~s/\.tsv$/.BiGG.bounds.tsv/;
}else{
	$OutputFile.=".BiGG.bounds.tsv";
}

###############################################################################

print STDERR "\n";
print STDERR "Input Exchange Reaction File: $InputFile\n";
print STDERR "Translation Map: $MapFile\n";
print STDERR "Output File: $OutputFile\n";
print STDERR "\n";

###############################################################################

my %map;
my %degenerate_hash;
my $num_map_entries=0;

open(MAP_FH, "<$MapFile") || die "Could not open $MapFile for reading\n";
while(<MAP_FH>){
	chomp;
	my ($key, $val)=split /\t/, $_;

	if($val ne ""){
		if(!defined($map{$key})){
			$map{$key}=$val;
		}else{
			print STDERR "Duplicate key/val mapping: $key -> $val (Also mapped to: $map{$key})\n";
			$degenerate_hash{$key}=1;
		}
	}else{
		$map{$key}="UNMAPPED";
	}
	$num_map_entries++;
}
close(MAP_FH);

print STDERR "Num Map Entries Read: $num_map_entries\n";
print STDERR "\n";

###############################################################################

delete($map{""});

###############################################################################

sub parse_exchange_equation{
	my $equation=shift;

	my ($comp_str, $equation_str)=split " : ", $equation;
	
	if($equation_str eq ""){
		$equation_str = $comp_str;
	}

	my ($met_id)=split " <==>", $equation_str;

	return($met_id);
}

open(RXN_FH, "<$InputFile") || die "Could not open $InputFile for reading.\n";
open(OUT_FH, ">$OutputFile") || die "Could not open $OutputFile for writing.\n";

while(<RXN_FH>){
	chomp;
	
	my @fields=split "\t", $_;

	if($fields[0] eq "abbreviation"){
		next;
	}else{
		my $equation=$fields[2];
		my $lowbnd=$fields[5];
		my $uppbnd=$fields[6];
		my $description=$fields[1];

		my $metabolite=parse_exchange_equation($equation);

		#print STDERR "$equation / '$metabolite'\n";

		my $bigg_met_id=$map{$metabolite};
		if($bigg_met_id eq ""){
			$bigg_met_id=$metabolite;
		}

		my $line=join "\t", ($bigg_met_id, $lowbnd, $uppbnd, "$description / $metabolite");
		print OUT_FH "$line\n";
	}

}

close(RXN_FH);


###############################################################################

print STDERR "Done.\n";
