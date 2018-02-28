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
		-r [reverse flag, by default best is small, else reversing causses best to be the largest ]

	Reads in input file and based on the column you specify, picks all the records with a best value.
	Output will be written to <input file>.best.  The id (which groups what the best should
	be chosen from), should be column 0.

	If there is more than one equally best, then among the best, one is randomly selected.
	This is fine for large datasets where you are just trying to a get distribution of
	where hits were, but this results should not be considered definitive.

";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_c;
my $large_is_best=defined($opt_r);

my $ext;
if(defined($large_is_best)){
	$ext="best_large";
}else{
	$ext="best_small";
}

###############################################################################

open(IN_FH, "<$file") || die "Could not open $file\n";
open(OUT_FH, ">$file\.$ext") || die "Could not open $file\.$ext";

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

sub pick_best{
	my $sorted_arr_ref=shift;
	my $num_entries=$#{$sorted_arr_ref}+1;
	
	# The first entry is the best, becasue it's already sorted
	my $best_entry=${$sorted_arr_ref}[0];
	my @best_entry_fields=split /\t/, $best_entry;
	my $best_val=$best_entry_fields[$column]+0;
	#print "$best_val\n";

	# Figure out how many in the top are tied
	my $i;
	for($i=1; $i<$num_entries; $i++){
		my $cur_entry=${$sorted_arr_ref}[$i];
		my @cur_entry_fields=split /\t/, $cur_entry;
		my $cur_val=$cur_entry_fields[$column]+0;

		if($cur_val!=$best_val){
			last;
		}
	}

	# Randomly pick a number between 0 and (Num Tied - 1)
	my $num_best=$i;
	my $rand_best_idx=int(rand()*$num_best);

	# Return randomly selected best
	return(${$sorted_arr_ref}[$rand_best_idx]);
}

sub dump_best{
	my $arr_ref=shift;

	# If no entries, return nothing
	if($#{$arr_ref}==-1){
		return;
	}

	# Sort entries
	my @sorted_arr;
	if($large_is_best){
		@sorted_arr=sort large_to_small @{$arr_ref};
	}else{
		@sorted_arr=sort small_to_large @{$arr_ref};
	}

	# Randomly pick a best, if there is a tie
	my $rand_best_entry=pick_best(\@sorted_arr);

	#for(my $i=0; $i<=$#{$arr_ref}; $i++){
	#	print ${$arr_ref}[$i] . "\n";
	#}
	#print "\n";
	#print $rand_best_entry . "\n";
	#print "\n\n\n";

	print OUT_FH "$rand_best_entry\n";

}

