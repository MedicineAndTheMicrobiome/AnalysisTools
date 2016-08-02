#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_t);

getopts("t:");

my $usage = "
	usage:
	$0
		-t <EC dat file>

	Reads in EC enzyme.dat file and generates a simplified
	table of ID and description.

	Output goes to STDOUT.
	
";

if(!defined($opt_t)){
	die $usage;
}

my $ec_dat_filename=$opt_t;

###############################################################################

my $cur_AC="";
my $cur_DE="";

open(FH, "<$ec_dat_filename") || die "Could not open $ec_dat_filename\n";

my $i=0;
while(<FH>){
	chomp $_;
	
	if($_=~/^ID   /){

		# Extract EC ID
		my ($tag, $rec)=split "   ", $_;
		$cur_AC=$rec;

	}elsif($_=~/^DE   /){

		# Extract description	
		my ($tag, $rec)=split "   ", $_;
		$rec=~s/\.$//;
		$cur_DE=$rec;

	}elsif($_=~/^\/\//){

		# Dumping record out.

		my $outstr=join "\t", ($cur_AC, $cur_DE);

		if($outstr ne "\t"){
			print "$outstr\n";
		}

		$cur_AC=$cur_DE="";

		# Heart beat
		if(!($i%1000)){
			print STDERR "$i Records Processed.\n";
		}
		$i++;
	}
}
close(FH);

###############################################################################

print STDERR "Done.\n";
