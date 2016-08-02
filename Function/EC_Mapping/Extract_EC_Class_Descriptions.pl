#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_t);

getopts("t:");

my $usage = "
	usage:
	$0
		-t <EC class txt file>

	Extract the enzyme classes.

	Output goes to STDOUT.
	
";

if(!defined($opt_t)){
	die $usage;
}

my $ec_dat_filename=$opt_t;

###############################################################################

my $L1="";
my $L2="";
my $L3="";

open(FH, "<$ec_dat_filename") || die "Could not open $ec_dat_filename\n";

my $i=0;
while(<FH>){
	chomp $_;
	
	if($_=~/ *(\d*)\. *(\d*-*)\. *(\d*-*)\. *(\d*-*) +(.+)\./){
		#print "$1.$2.$3.$4\t$5\n";

		if($1 ne "-" && $2 eq "-" && $3 eq "-"){
			$L1=$5;
			$L2="";
			$L3="";
		}
		if($1 ne "-" && $2 ne "-" && $3 eq "-"){
			$L2=$5;
			$L3="";
		}
		if($1 ne "-" && $2 ne "-" && $3 ne "-"){
			$L3=$5;
		}

		my $outstr="$1.$2.$3.$4\t$L1 :: $L2 :: $L3\n";

		$outstr=~s/ :: $//;
		$outstr=~s/ :: $//;
		
		print $outstr;
	}

}
close(FH);

###############################################################################

print STDERR "Done.\n";
