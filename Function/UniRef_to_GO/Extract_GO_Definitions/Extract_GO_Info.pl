#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_m $opt_o);

getopts("m:o:");

my $usage = "
	usage:
	$0
		-m <GO OBO File>
		-o <Output file name>

	This script will extract GO Definitions into a tab-delimited file.

";

if(!defined($opt_m) || !defined($opt_o)){ 
	die $usage;
}

my $dat_file=$opt_m;
my $output_fn=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "GO OBO File: $dat_file\n";
print STDERR "Output File: $output_fn\n";
print STDERR "\n";

###############################################################################

open(IN_FH, "<$dat_file") || die "Could not open dat file $dat_file\n";
open(OUT_FH, ">$output_fn") || die "Could not open output file $output_fn\n";

my $NEW_REC_KEYWORD="[Term]";
my $MULTILINE=0;

my $id;
my $name;
my @is_a; 

while(<IN_FH>){
	chomp;

	if(($_ eq $NEW_REC_KEYWORD) && ($id ne "")){
	
		my $num_isas=$#is_a+1;
		my $is_a_string;

		if($num_isas==0){
			$is_a_string="GO:-1";
			print OUT_FH "$id\t$is_a_string\t$name\n";
		}else{
			if(!$MULTILINE){
				$is_a_string = join ";", @is_a;
				print OUT_FH "$id\t$is_a_string\t$name\n";
			}else{
				for(my $i=0; $i<$num_isas; $i++){
					print OUT_FH "$id\t$is_a[$i]\t$name\n";
				}
			}
		}

		$id="";
		$name="";
		@is_a=();		

	}else{

		if($_=~/([a-zA-Z0-9_]+): (.+)/){
			my $key=$1;
			my $value=$2;
			if($key eq "id"){
				$id=$value;
			}elsif($key eq "name"){
				$name=$value;
			}elsif($key eq "is_a"){
				my ($id, $rel_str)=split / ! /, $value;
				push @is_a, $id;
			}
		
		}
	}

}

close(IN_FH);
close(OUT_FH);

###############################################################################

print STDERR "Done.\n";
