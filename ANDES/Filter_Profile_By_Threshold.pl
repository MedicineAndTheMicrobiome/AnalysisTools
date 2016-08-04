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

	One of the following:
	-t <base count threshold>
	-p <percentage threshold, ie. 5, is 5%> 
	-b <confidence of presence, ie. 95 is 95% sure>

	Reads in a profile file and sets the base count to zero if its
	value does not equal or exceed the threshold set.  It's essentially
	a noise filter.  Once the counts have been reset, if a major allele
	can not be detected, then an N will be used.
	
";

if(
	!defined($opt_i) ||
	!defined($opt_o) ||
	(!defined($opt_t) && !defined($opt_p) && !defined($opt_b)) ||
	(defined($opt_t) + defined($opt_p) + defined($opt_b))>1
){
	die $usage;
}

my $input_prof_fname=$opt_i;
my $output_prof_fname=$opt_o;
my $base_count_thres=$opt_t;
my $percentage_thres=$opt_p;
my $binomial_confidence=$opt_b;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

my $filtered_prof_arr_ref;

if(defined($base_count_thres)){
	$filtered_prof_arr_ref=PositionProfiles::filter_prof_basecount($prof_arr_ref, $base_count_thres);
}elsif(defined($percentage_thres)){
	$filtered_prof_arr_ref=PositionProfiles::filter_prof_percentage($prof_arr_ref, $percentage_thres/100.0);
}elsif(defined($binomial_confidence)){
	$filtered_prof_arr_ref=PositionProfiles::filter_prof_binomial($prof_arr_ref, $binomial_confidence/100.0);
}

# Remove positions where there is nothing but gaps
$filtered_prof_arr_ref=PositionProfiles::remove_pos_with_all_gaps($filtered_prof_arr_ref);

my $major_nucs=PositionProfiles::compute_mode_residue_for_profile($filtered_prof_arr_ref);
PositionProfiles::output_prof($output_prof_fname, $filtered_prof_arr_ref, $major_nucs);

