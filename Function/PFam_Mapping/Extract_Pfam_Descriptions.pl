#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_t);

getopts("t:");

my $usage = "
	usage:
	$0
		-t <Pfam mapping file>

	Reads in pdb_pfam_mapping.txt file

	Output goes to STDOUT.
	
";

if(!defined($opt_t)){
	die $usage;
}

my $pfam_dat_filename=$opt_t;

###############################################################################


open(FH, "<$pfam_dat_filename") || die "Could not open $pfam_dat_filename\n";

# Skip over header line
$_=<FH>;

my $i=0;
my $ACC_COL=4;
my $DESC_COL=6;
while(<FH>){
	chomp $_;

	my @val=split "\t", $_;

	my ($clean_acc, $version)=split /\./, $val[$ACC_COL];
	my $desc=$val[$DESC_COL];

	print "$clean_acc\t$desc\n";
}
close(FH);

###############################################################################

print STDERR "Done.\n";
