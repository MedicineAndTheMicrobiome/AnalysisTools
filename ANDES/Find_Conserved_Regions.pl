#!/usr/bin/env perl

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_p $opt_l $opt_i $opt_o);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("p:l:i:o:");
my $usage = "usage: 
$0 
	-p <input profile>
	-l <region/window length, e.g. 40 (bp)>
	-i <min identity, e.g. 95 (%)>
	-o <output root>

This script will read in the specified input profile (-p) and
look for regions of length (-l) that meet or exceed the specified
percent idenity cutoff (-i).  

The output will be:
	1.) Contiguous FASTA file with a record for each region that was found to exceed the -l length
	2.) Overlapping FASTA file, where each record is an overlapping sequenced windowed from the contiguous region
	3.) Begin/end/sequence in 0-space-based coordinates of the contiguous regions

";

if(
	!defined($opt_p) ||
	!defined($opt_l) ||
	!defined($opt_i) ||
	!defined($opt_o)
){
	die $usage;
}

my $in_prof=$opt_p;
my $regions_len=$opt_l;
my $min_ident=$opt_i;
my $output_root=$opt_o;

print STDERR "Profile: $in_prof\n";
print STDERR "Regions Length: $regions_len\n";
print STDERR "Min. Identity: $min_ident\n";
print STDERR "Output Root: $output_root\n";

###############################################################################

sub find_regions{
	my $mask_arr_ref=shift;
	my $mask_len=$#{$mask_arr_ref}+1;
	my @mask=@{$mask_arr_ref};

	my @positions_arr;
	my @lengths_arr;

	foreach my $pos(@mask){
		print STDERR $pos;
	}
	print STDERR "\n";

	for(my $i=0; $i<$mask_len; $i++){
		if($mask[$i]==1){

			# Store position of region
			push @positions_arr, $i;

			# Compute length of region
			my $length=0;
			while($mask[$i]==1 && $i<$mask_len){
				$length++;
				$i++;
			}

			# Store length of region
			push @lengths_arr, $length;
		}
	}

	return(\@positions_arr, \@lengths_arr);	
}

###############################################################################

#my @test=(1,0,0,0,1,1,1,0,1,1,1,1,1);
#my ($p,$r)=find_regions(\@test);
#for(my $i=0; $i<=$#{$r}; $i++){
#	print "${$p}[$i]\t${$r}[$i]\n";
#}
#exit;

my ($major_nuc_ref, $prof_ref)=PositionProfiles::load_prof($in_prof);
my $prof_length= PositionProfiles::profile_arr_length($prof_ref);
printf ("Input Profile Length: %i\n", $prof_length);

my $mininum_required_depth=100/(100-$min_ident);
print STDERR "Minimum Coverage Required: $mininum_required_depth bps\n";


my ($raw_null, $raw_valid)=PositionProfiles::count_null_positions($prof_ref);
print STDERR "Input: Null: $raw_null  Valid: $raw_valid\n";

my $zero_masked_prof=PositionProfiles::filter_prof_depth($prof_ref, $mininum_required_depth);
#PositionProfiles::print_prof($zero_masked_prof);
printf("Depth Filtered: %i\n", PositionProfiles::profile_arr_length($prof_ref));

my ($dfilt_null, $dfilt_valid)=PositionProfiles::count_null_positions($zero_masked_prof);
print STDERR "Depth Filtered: Null: $dfilt_null  Valid: $dfilt_valid\n";



my $norm_prof=PositionProfiles::normalize_profile($zero_masked_prof);
#PositionProfiles::print_prof($norm_prof);
printf("Normalized Length: %i\n", PositionProfiles::profile_arr_length($norm_prof));

my $gap_free_norm_prof=PositionProfiles::filter_gaps($norm_prof);
printf("Gap Filtered Length: %i\n", PositionProfiles::profile_arr_length($gap_free_norm_prof));

my $mask_arr_ref=PositionProfiles::mask_positions_over_threshold($gap_free_norm_prof, $min_ident/100);
my $gap_free_major_seq=PositionProfiles::compute_mode_residue_for_profile($gap_free_norm_prof);

my ($pos_ref, $len_ref)=find_regions($mask_arr_ref);

my $num_regions=$#{$pos_ref}+1;

my @filtered_pos;
my @filtered_len;

for(my $i=0; $i<$num_regions; $i++){
	if(${$len_ref}[$i] > $regions_len){
		push @filtered_len, ${$len_ref}[$i];
		push @filtered_pos, ${$pos_ref}[$i];
	}
}

output_mask("$output_root.$min_ident", $mask_arr_ref);

output_regions("$output_root.$min_ident", \@filtered_pos, \@filtered_len, $gap_free_major_seq);

output_contiguous_fasta("$output_root.$min_ident.ctg", \@filtered_pos, \@filtered_len, $gap_free_major_seq);

output_overlapping_fasta("$output_root.$min_ident.ovl", \@filtered_pos, \@filtered_len, $gap_free_major_seq, $regions_len);


###############################################################################

sub output_mask{
	my $fname=shift;
	my $mask_ref=shift;
	my $len=$#{$mask_ref}+1;

	open(FH, ">$fname\.mask") || die "Could not open $fname\.mask for writing\n";

	for(my $i=0; $i<$len; $i++){
		print  FH "$i\t${$mask_ref}[$i]\n";
	}	

	close(FH);
}

###############################################################################

sub output_contiguous_fasta{
	my $fname=shift;
	my $pos_ref=shift;
	my $len_ref=shift;
	my $seq=shift;

	my @glob_seq_arr=split //, $seq;

	print STDERR "Outputing FASTA file...\n";

	open(FH, ">$fname\.fasta") || die "Could not open $fname\.fasta for writing.\n";

	for(my $i=0; $i<=$#{$pos_ref}; $i++){
		my $begin=${$pos_ref}[$i];
		my $length= ${$len_ref}[$i];
		my $end=$begin+$length;


		my @reg_seq_arr=@glob_seq_arr[ $begin .. $begin+$length-1];
		my $reg_seq=join "", @reg_seq_arr;
		print FH ">$fname\.$begin-$end\.$length\n$reg_seq\n";
	}

	close(FH);


}

###############################################################################


sub output_overlapping_fasta{
	my $fname=shift;
	my $pos_ref=shift;
	my $len_ref=shift;
	my $seq=shift;
	my $regions_len=shift;

	print STDERR "Outputing FASTA file...\n";

	open(FH, ">$fname\.fasta") || die "Could not open $fname\.fasta for writing.\n";

	my @glob_seq_arr=split //, $seq;

	for(my $i=0; $i<=$#{$pos_ref}; $i++){
		my $begin=${$pos_ref}[$i];
		my $length= ${$len_ref}[$i];

		for(my $j=0; $j<=($length-$regions_len); $j++){
			my $window_begin=$begin+$j;
			my $window_end=$window_begin+$regions_len;
			my @reg_seq_arr=@glob_seq_arr[ $window_begin .. $window_begin+$regions_len-1];
			my $reg_seq=join "", @reg_seq_arr;
			print FH ">$fname\.$window_begin-$window_end\.$length\n$reg_seq\n";
		}
	}

	close(FH);


}

###############################################################################

sub output_regions{
	my $fname=shift;
	my $pos_ref=shift;
	my $len_ref=shift;
	my $seq=shift;

	print STDERR "Outputing Regions file...\n";
	open(FH, ">$fname\.regions") || die "Could not open $fname\.regions for writing.\n";

print "$seq\n";
	my @tmp_seq=split //, $seq;

	for(my $i=0; $i<=$#{$pos_ref}; $i++){
		my @reg_seq_arr=@tmp_seq[ ${$pos_ref}[$i] ..  ${$pos_ref}[$i]+${$len_ref}[$i]-1];
		my $reg_seq=join "", @reg_seq_arr;
		print FH "${$pos_ref}[$i]\t${$len_ref}[$i]\t$reg_seq\n";
	}

	close(FH);


}

###############################################################################

