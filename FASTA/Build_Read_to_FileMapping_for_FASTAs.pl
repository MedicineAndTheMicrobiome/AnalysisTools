#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_o $opt_s $opt_e);
use File::Basename;

getopts("o:se");
my $usage = "usage: 
$0 
	-o <output root>
	[-s (strip path option)]
	[-e (strip fasta extension option)]
	<list of fasta files>...

	This script will read in all the fasta file and produce a mapping
	of read ids to output file names.

	<read_id>\\t<filename>\\n
	<read_id>\\t<filename>\\n
	<read_id>\\t<filename>\\n
	...
	<read_id>\\t<filename>\\n

	By default, the entire path and filename will be used as the <filename>
	Use the -s or -e to strip the directory path or .fa/.fasta extension,
	respectively.

";

if(!(
	defined($opt_o))){
	die $usage;
}

###############################################################################

my $output=$opt_o;
my $strip_path=defined($opt_s);
my $strip_ext=defined($opt_e);

print STDERR "Output root: $output\n";
open(OUT, ">$output") || die "Could not open $output\n";

while(my $filename=shift){
	print STDERR "Working on: $filename\n";

	my $src=$filename;

	if($strip_ext){
		$src=~s/\.fasta$//;
		$src=~s/\.fa$//;
	}
	if($strip_path){
		($src)=fileparse($src);
	}

	open(FH, "<$filename") || die "Could not open $filename\n";
	
	while(<FH>){
		if($_=~/^>(\S+)/){
			chomp;
			print OUT "$1\t$src\n";
		}
	}

	close(FH);
}

#------------------------------------------------------------------------------

