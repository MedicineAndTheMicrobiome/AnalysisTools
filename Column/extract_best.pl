#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c $opt_r $opt_t);

getopts("i:c:t:r");

my $usage = "
	usage:
	$0
		-i <input file>
		-c <column>
		[-t <tie breaker column>]
		-r [reverse flag, best is small, else best is largest column value]

	Reads in input file and based on the column you specify, picks the one with the best value.
	Output will be written to <input file>.best.  The id (which groups what the best should
	be chosen from), should be column 0.

	You can use the -t option to make sure that the best is selected consistently if there is more
	than one.  The tie breaker is always in ascending order.  Make sure the tie breaker is
	unique for a given set of groupings.

";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_c;
my $small_is_best=$opt_r;
my $tie_breaker_column=$opt_t;

if(!defined($opt_t)){
	$tie_breaker_column=0;
}

###############################################################################

open(IN_FH, "<$file") || die "Could not open $file\n";
open(OUT_FH, ">$file\.best") || die "Could not open $file\.best";

my $cur_id;
my $pre_id;
my @array;
my @acc_arr;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;
	$cur_id=$array[0];
	if($cur_id eq $pre_id){
		push @acc_arr, $_;
	}else{
		dump_best(\@acc_arr);
		@acc_arr=();
		push @acc_arr, $_;
	}
	$pre_id=$cur_id;
}
dump_best(\@acc_arr);

close(IN_FH);
close(OUT_FH);

###############################################################################

sub small_to_large{
	my @a_arr=split /\t/, $a;
	my @b_arr=split /\t/, $b;
	
	($a_arr[$column]+0)<=>($b_arr[$column]+0)
	or
	($a_arr[$tie_breaker_column]+0)<=>($b_arr[$tie_breaker_column]+0);
}

sub large_to_small{
	my @a_arr=split /\t/, $a;
	my @b_arr=split /\t/, $b;
	
	($b_arr[$column]+0)<=>($a_arr[$column]+0)
	or
	($a_arr[$tie_breaker_column]+0)<=>($b_arr[$tie_breaker_column]+0);
}

sub dump_best{
	my $arr_ref=shift;
	if($#{$arr_ref}==-1){
		return;
	}
	my @sorted_arr;

	if($small_is_best){
		@sorted_arr=sort small_to_large @{$arr_ref};
	}else{
		@sorted_arr=sort large_to_small @{$arr_ref};
	}
	print OUT_FH "$sorted_arr[0]\n";
}

