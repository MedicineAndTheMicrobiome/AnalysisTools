#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i);

getopts("i:");

my $usage = "
	usage:
	$0
		-i <Input Reaction Filename>

	This script wil read in a reactions file and report which
	exchange reactions are:
	
		1.) Fully open
		2.) In only
		3.) Out only
		4.) Shutdown
		5.) Tweaked

";

if(!defined($opt_i)){
	die $usage;
}

my $InputFile=$opt_i;
my $outputfile_root=$InputFile;
$outputfile_root=~s/_react\.tsv$//;
my $outputfile="$outputfile_root.exch_flux_restr.tsv";

###############################################################################

print STDERR "Reaction Input: $InputFile\n";
print STDERR "Output File: $outputfile\n";

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

my $fh;
my $temp_filename;

my $ID_COLUMN=0;
my $DESC_COLUMN=1;
my $REACTION_COLUMN=2;
my $LB_COLUMN=5;
my $UB_COLUMN=6;

my %ex_full;
my %ex_in;
my %ex_out;
my %ex_shutdown;
my %ex_tweaked;

my $MAX=1000;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;

	$#array=9;
	my $outstr=join "\t", @array;

	# Grab reaction column
	my $reaction=$array[$REACTION_COLUMN];

	if($_=~/<==>/){
		my ($lhs, $rhs)=split /<==>/, $reaction;
		$rhs=~s/\s+//g;

		if($rhs eq ""){
			# print "$outstr\n";

			my $desc=$array[$DESC_COLUMN];
			my $lb=$array[$LB_COLUMN];
			my $ub=$array[$UB_COLUMN];

			if($lb==0 && $ub==0){
				$ex_shutdown{$desc}=$reaction;
			}elsif($lb==0 && $ub==$MAX){
				$ex_out{$desc}=$reaction;
			}elsif($lb==-$MAX && $ub==0){
				$ex_in{$desc}=$reaction;
			}elsif($lb==-$MAX && $ub==$MAX){
				$ex_full{$desc}=$reaction;
			}else{
				$ex_tweaked{$desc}="$reaction\t$lb\t$ub";
			}

		}
	}

	# Output pulse
	$input_lines++;
}

close(IN_FH);

###############################################################################

open(OUT_FH, ">$outputfile") || die "Could not open $outputfile for writing.\n";

print OUT_FH "\nExchanges Shutdown (No Flux Allowed):\n";
if(scalar keys %ex_shutdown){
	foreach my $reaction(sort keys %ex_shutdown){
		print OUT_FH "\t$reaction:\t$ex_shutdown{$reaction}\n";
	}
}else{
	print OUT_FH "\t<none>\n";
}

print OUT_FH "\nExchanges Full (Unrestricted Flux):\n";
if(scalar keys %ex_full){
	foreach my $reaction(sort keys %ex_full){
		print OUT_FH "\t$reaction:\t$ex_full{$reaction}\n";
	}
}else{
	print OUT_FH "\t<none>\n";
}

print OUT_FH "\nExchanges In Only (Only Consumption Allowed):\n";
if(scalar keys %ex_in){
	foreach my $reaction(sort keys %ex_in){
		print OUT_FH "\t$reaction:\t$ex_in{$reaction}\n";
	}
}else{
	print OUT_FH "\t<none>\n";
}

print OUT_FH "\nExchanges Out Only (Only Production Allowed):\n";
if(scalar keys %ex_out){
	foreach my $reaction(sort keys %ex_out){
		print OUT_FH "\t$reaction:\t$ex_out{$reaction}\n";
	}
}else{
	print OUT_FH "\t<none>\n";
}

print OUT_FH "\nExchanges Tweaked (Non-Standard Bounds Defined):\n";
if(scalar keys %ex_tweaked){
	foreach my $reaction(sort keys %ex_tweaked){
		print OUT_FH "\t$reaction:\t$ex_tweaked{$reaction}\n";
	}
}else{
	print OUT_FH "\t<none>\n";
}


print OUT_FH "\n";
close(OUT_FH);

###############################################################################

print STDERR "\nDone.\n";
