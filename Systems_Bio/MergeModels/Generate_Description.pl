#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Basename;
use vars qw ($opt_r $opt_m $opt_o);

getopts("r:m:o:");

my $usage = "
	usage:
	$0
		-r <reaction.tsv>
		-m <met.tsv>
		-o <output directory>
";

if(!defined($opt_r) || 
	!defined($opt_m)
){
	die $usage;
}


my $reactname=$opt_r;
my $metname=$opt_m;
my $output_dir=$opt_o;

if(!defined($output_dir)){
	$output_dir=".";
}

print STDERR "Reaction Filename: $reactname\n";
print STDERR "Metabolites Filename: $metname\n";
print STDERR "Output Directory: $output_dir\n";
print STDERR "\n";

###############################################################################
# Get root for description filename

my ($name, $path, $suffix)=fileparse($reactname);

my $descname=$name;
$descname=~s/_react\.tsv$//;
my $model_id=$descname;
$descname.="_desc.tsv";

print STDERR "Description Filename: $descname\n";

###############################################################################
# Get num metabolites

my $lc=`wc -l $metname`;
$lc=~s/^\s+//g;
my $nmet;
if($lc=~/^(\d+) /){
	$nmet=$1-1;	
}else{
	die "Error:  Problem getting line count from wc -l $metname\n";
}
print STDERR "Num metabolites found: $nmet\n";

###############################################################################
# Get num reactions

my $lc=`wc -l $reactname`;
$lc=~s/^\s+//g;
my $nreact;
if($lc=~/^(\d+) /){
	$nreact=$1-1;	
}else{
	die "Error:  Problem getting line count from wc -l $reactname\n";
}
print STDERR "Num reactions found: $nreact\n";

###############################################################################
# Get Compartments

open(FH, "<$reactname") || die "Could not open $reactname\n";
<FH>; #Ignore header
my $COMPART_COL=4;
my %comp_hash;
while(<FH>){
	my @fields=split "\t", $_;
	my @comparts=split ",", $fields[$COMPART_COL];
	foreach my $comp (@comparts){
		$comp=~s/^\s+//g;
		$comp=~s/\s+$//g;
		$comp_hash{$comp}=1;	
	}
}
close(FH);

my @foundcomp=sort keys %comp_hash;
my $compartments=join ",", @foundcomp;

for(my $i=0; $i<=$#foundcomp; $i++){
	$foundcomp[$i]="[" . $foundcomp[$i] . "]";
}
my $abbreviations=join ",", @foundcomp;

###############################################################################

my @headers=(
	"name",
	"id",
	"description",
	"compartment",
	"abbreviation"
#	"Nmetabolites",
#	"Nreactions",
#	"Ngenes",
#	"Nnnz"
);

my @values=(
	$model_id,
	$model_id,
	$model_id,
	$compartments,
	$abbreviations
#	$nmet,
#	$nreact,
#	"",
#	""	
);

###############################################################################

open(OUT_FH, ">$output_dir/$descname") || die "Could not open $$output_dir/$descname\n";

my $outline=join "\t", @headers;
print OUT_FH "$outline\n";

$outline=join "\t", @values;
print OUT_FH "$outline\n";

close(OUT_FH);

###############################################################################

print STDERR "\n";
print STDERR "Done.\n";
