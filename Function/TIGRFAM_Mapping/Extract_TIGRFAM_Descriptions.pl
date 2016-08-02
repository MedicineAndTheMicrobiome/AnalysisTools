#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_t);

getopts("t:");

my $usage = "
	usage:
	$0
		-t <TIGRFAM info file>

	Reads in TIGRFAM info file (that was created by concatenation)
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
	
	if($_=~/^AC  /){

		# Extract EC ID
		my ($tag, $rec)=split "  ", $_;
		$cur_AC=$rec;

	}elsif($_=~/^DE  /){

		# Extract description	
		my ($tag, $rec)=split "  ", $_;
		$cur_DE=$rec;

		if($cur_AC ne ""){
			print "$cur_AC\t$cur_DE\n";
		}

	}
}

close(FH);

###############################################################################

print STDERR "Done.\n";
