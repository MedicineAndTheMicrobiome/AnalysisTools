#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c $opt_r);

getopts("i:c:r");

my $usage = "
	usage:
	$0
		-i <input file>
		-c <column>
		-r [reverse flag, best is small, else best is largest column value]

	Reads in input file and based on the column you specify, picks all the records with a best value.
	Output will be written to <input file>.all_best.  The id (which groups what the best should
	be chosen from), should be column 0.

";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_c;
my $small_is_best=$opt_r;

###############################################################################

open(IN_FH, "<$file") || die "Could not open $file\n";
open(OUT_FH, ">$file\.all_best") || die "Could not open $file\.all_best";

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
	
	($a_arr[$column]+0)<=>($b_arr[$column]+0);
}

sub large_to_small{
	my @a_arr=split /\t/, $a;
	my @b_arr=split /\t/, $b;
	
	($b_arr[$column]+0)<=>($a_arr[$column]+0);
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

	my $more_best=1;
	my $i=0;
	my $best_val=(split /\t/, $sorted_arr[0])[$column];
	while($more_best){
		my @values=split /\t/, $sorted_arr[$i];
		if(@values[$column]==$best_val){
			print OUT_FH "$sorted_arr[$i]\n";
		}else{
			$more_best=0;
		}
		$i++;
	}
}

