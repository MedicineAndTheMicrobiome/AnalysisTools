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
use vars qw($opt_p $opt_w $opt_s);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

my $WINDOW=10;
my $STEP=10;

getopts("p:w:s:");
my $usage = "usage: 
$0 

	-p <profile>
	[-w <GC window size, $WINDOW bases>]
	[-s <step size, $STEP bases>]

	Reads in a profile, and generates data for a scatter plot.
	The first column (x) is the %%GC content, for a window.
	The second column (y) is the amount of Shannon Normalized Entropy (Evenness), between 0 and 1, experienced in that window.
	
	0 variation means only one base is present at a position.
	1 variation means all bases are represented equally at a position.

";

if(!defined($opt_p)){
	die $usage;
}

my $profile_name=$opt_p;
my $window_size=$WINDOW;

if(defined($opt_w)){
	$window_size=$opt_w;
}

my $step_size=$STEP;
if(defined($opt_s)){
	$step_size=$opt_s;
}

###############################################################################

my ($seq_ref, $prof_ref)=PositionProfiles::load_prof($profile_name);
my $normalized_prof_ref=PositionProfiles::normalize_profile($prof_ref);

#PositionProfiles::print_prof($normalized_prof_ref);

my ($variation_arr_ref)=generate_positional_variance($normalized_prof_ref);
my ($gc_content_ref, $variation_ref)=generate_gc_scatter_pairs($window_size,$step_size,$seq_ref,$variation_arr_ref);

my $scatter_plot_count=$#{$gc_content_ref}+1;
for(my $i=0; $i<$scatter_plot_count; $i++){
	print "${$gc_content_ref}[$i]\t${$variation_ref}[$i]\n";
}

###############################################################################

sub generate_gc_scatter_pairs{
	my $window_size=shift;
	my $step_size=shift;
	my $sequence_arr_ref=shift;
	my $variation_arr_ref=shift;

	# Empty arrays to keep track of window samples 
	my @gc_arr;
	my @var_arr;

	# Step through each window
	my $length=$#{$variation_arr_ref}+1;
	for(my $i=0; $i<$length; $i+=$step_size){

		# Compute window extents
		my $samp_begin=$i;
		my $samp_end=$samp_begin+$window_size;

		# Make sure we don't have any out of bound errors
		if($samp_end>$length){
			last;
		}else{
			my $sum_gc=0;
			my $sum_var=0;
		
			# Step through positions in current window
			for(my $j=$samp_begin; $j<$samp_end; $j++){

				# Compute GC
				my $nuc=${$sequence_arr_ref}[$j];
				if($nuc eq "G" || $nuc eq "C"){
					$sum_gc++;
				}

				# Sum up variances
				$sum_var+=${$variation_arr_ref}[$j];
			}

			# Keep track of averages
			push @gc_arr, ($sum_gc/$window_size);
			push @var_arr, ($sum_var/$window_size);

		}
	}

	# Return parallel array of windowed sample results
	return(\@gc_arr, \@var_arr);
}

###############################################################################

sub generate_positional_variance{
	my $prof_ref=shift;
	my @variation;

	my $num_pos=$#{$prof_ref}+1;
	for(my $i=0; $i<$num_pos; $i++){
		$variation[$i]=PositionProfiles::normalized_shannon_entropy(${$prof_ref}[$i]);
	}

	return(\@variation);
}
