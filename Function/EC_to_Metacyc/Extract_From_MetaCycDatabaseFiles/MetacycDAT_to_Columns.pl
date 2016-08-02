#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_m $opt_f $opt_o);

getopts("m:f:o:");

my $usage = "
	usage:
	$0
		-m <MetaCyc DAT file>
		-f <Fields to extract, comma-separated>
		-o <Output file name>
";

if(!defined($opt_m) || !defined($opt_f) || !defined($opt_o)){ 
	die $usage;
}

my $dat_file=$opt_m;
my $fields=$opt_f;
my $output_fn=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "MetaCyc DAT File: $dat_file\n";
print STDERR "Fields: $fields\n";
print STDERR "Output File: $output_fn\n";

###############################################################################

my @fields=split ",", uc($fields);

my $num_fields=$#fields+1;
my %fields_hash;

print STDERR "\nTargeted Fields:\n";
foreach my $targeted_field(@fields){
	print STDERR "\t$targeted_field\n";
	$fields_hash{$targeted_field}=1;
}

my $UNIQUE_ID="UNIQUE-ID";

my %record_hash;
$fields_hash{$UNIQUE_ID}=1;
foreach my $field(@fields){
	@{$record_hash{$field}}=();
}


###############################################################################

open(IN_FH, "<$dat_file") || die "Could not open dat file $dat_file\n";
open(OUT_FH, ">$output_fn") || die "Could not open output file $output_fn\n";

while(<IN_FH>){
	chomp;
	if($_[0] eq "#"){
		next;
	}elsif($_ eq "//"){

		my @out_arr;
		foreach my $field(@fields){
			my $str=join ";", @{$record_hash{$field}};
			#print "$field: $str\n";
			push @out_arr, $str;
			#print "$field: $record_hash{$field}\n";
		}
		print OUT_FH (join "\t", @out_arr) . "\n";

		# New record
		# Clear variables
		foreach my $field(@fields){
			$#{$record_hash{$field}}=-1;
		}

	}else{
		my @line_arr=split " - ", $_;
		my $key=shift @line_arr;
		
		if(defined($fields_hash{$key})){
			my $value=join " - ", @line_arr;
			push @{$record_hash{$key}}, $value;
		}

	}

}

close(IN_FH);

###############################################################################

###############################################################################

print STDERR "Done.\n";
