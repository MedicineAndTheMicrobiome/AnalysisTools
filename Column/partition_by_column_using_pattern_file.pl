#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_i $opt_c $opt_p $opt_s);

getopts("i:c:p:s:");

my $MAX_CATEGORIES=30;

my $usage = "
	usage:
	$0
		-i <input file>
		-c <column, counting from 0>
		-p <file with list of patterns and group name>
		[-s <field separator>]

	Pattern file should look like:

	<group name>\\t<regular expression>\\n


	Reads in input file and based on the column you specify
	splits the file into multiple files using the value in that 
	column as a key.

	Output will be written to files named:

		<input file>.<group name>  
	
";

if(!defined($opt_i) || !defined($opt_c) || !defined($opt_p)){
	die $usage;
}

my $infile=$opt_i;
my $column=$opt_c;
my $pattern_file=$opt_p;
my $field_separator=$opt_s;

if(!defined($opt_s)){
	$field_separator="\t";
}

###############################################################################

my %patterns_hash;
my %groups_hash;

open(PATTERN_FH, "<$pattern_file") || die "Could not open pattern file\n";
while(<PATTERN_FH>){
	chomp;
	my ($group_id, $pattern)=split /\t/, $_;
	$patterns_hash{$pattern}=$group_id;
	$groups_hash{$group_id}=1;
}
close(PATTERN_FH);

###############################################################################

my @groups=keys %groups_hash;
my %filehandles;
for(my $i=0; $i<=$#groups; $i++){
	$filehandles{$groups[$i]}=new FileHandle ">$infile\.$groups[$i]";
}

###############################################################################

my @patterns=keys %patterns_hash;

###############################################################################

open(INFILE, "<$infile") || die "Could not open $infile\n";
while(<INFILE>){
	chomp;
	my @arr=split /$field_separator/, $_;
	for(my $i=0; $i<=$#patterns; $i++){
		my $curpat=$patterns[$i];
		if($arr[$column]=~/$curpat/){
			print {$filehandles{$patterns_hash{$curpat}}} "$_\n";
		}
	}
}
close(INFILE);


###############################################################################
