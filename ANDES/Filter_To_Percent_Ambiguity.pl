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
use vars qw($opt_i $opt_o $opt_t $opt_p $opt_b);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("i:o:t:p:b:");
my $usage = "usage: 
$0 
	-i <Input Profile>
	-o <Output Profile>
	-p <percentage ambiguity threshold, ie. 10, is 10%> 

	Reads in a profile file and then iteratively increases the percent cutoff
	for base representation until the overall percent ambiguity drops below
	the specified percent ambiguity threshold.
	
";

if(
	!(
	defined($opt_i) && defined($opt_o) && defined($opt_p)
	)
){
	die $usage;
}

my $input_prof_fname=$opt_i;
my $output_prof_fname=$opt_o;
my $ambiguity_threshold=$opt_p;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

my $threshold_increment=.1;
my $threshold=0;
my $filtered_prof_arr_ref;
my $perc_amb;

my @threshold_increment_arr=(10,1,.1,.01,.001,.0001,.00001);

my $threshold_increment;
my $used_threshold;
my $final_ambiguity;
foreach $threshold_increment(@threshold_increment_arr){

	do{

		# Apply threshold
		$threshold+=$threshold_increment;
		$used_threshold=$threshold;
		$filtered_prof_arr_ref=PositionProfiles::filter_prof_percentage($prof_arr_ref, $threshold/100.0);

		# Compute ambiguity
		my $amb_count=PositionProfiles::count_ambiguous_positions_across_profile($filtered_prof_arr_ref);
		my $prof_length=PositionProfiles::profile_arr_length($prof_arr_ref);	

		$perc_amb=100*$amb_count/$prof_length;
		$perc_amb=sprintf("%3.4f", $perc_amb);

		$final_ambiguity=$perc_amb;
		print STDERR "$amb_count of $prof_length are ambiguous ($perc_amb%) at $threshold% cutoff.\n";

	# Stop when too many ambiguities filtered
	}while($perc_amb>$ambiguity_threshold && $threshold<100);

	# Backoff, and increase granularity
	$threshold-=$threshold_increment;
}


# Remove positions where there is nothing but gaps
$filtered_prof_arr_ref=PositionProfiles::remove_pos_with_all_gaps($filtered_prof_arr_ref);

# Prepare and output filtered profile
my $major_nucs=PositionProfiles::compute_mode_residue_for_profile($filtered_prof_arr_ref);
PositionProfiles::output_prof($output_prof_fname, $filtered_prof_arr_ref, $major_nucs);

print STDERR "Done.  \nFinal threshold: $used_threshold%  Final percent ambiguity: $final_ambiguity%.\n";



