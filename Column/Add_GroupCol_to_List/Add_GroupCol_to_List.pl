#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_g);

getopts("g:");

my $usage = "
	usage:
	$0
		-g <group name>

	This script will read in a file through STDIN and append
	the <group name> to the next column.

	Example usage:

	cat var.lst | $0 -g snps > snp.grouping
	lch metadata.tsv | $0 -g snps > snp.grouping
	
";



if(!defined($opt_g) && (-t STDIN) && $ARGV[0] eq ""){
	die $usage;
}

my $Group=$opt_g;


###############################################################################

print STDERR "Running $0\n";

my $fh;

while(<STDIN>){

	chomp $_;
	my $inline=$_;
	print STDOUT "$inline\t$Group\n";

}

print STDERR "Done.\n";
