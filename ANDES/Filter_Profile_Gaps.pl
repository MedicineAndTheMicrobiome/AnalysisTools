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
use vars qw($opt_i $opt_o $opt_t);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("i:o:t:");
my $usage = "usage: 
$0 
	-i <Input Profile>
	-o <Output Profile>
	-t <max percentage threshold, ie. .75, is 75%> 

	Reads in the input file and then removes the position if 
	the percentage of gaps at that position is greater than
	the max percentage threshold.
	
";

if(
	!defined($opt_i) ||
	!defined($opt_o) ||
	!defined($opt_t) 
){
	die $usage;
}

my $input_prof_fname=$opt_i;
my $output_prof_fname=$opt_o;
my $percentage_thres=$opt_t;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

my $filtered_prof_arr_ref;
$filtered_prof_arr_ref=PositionProfiles::filter_gaps($prof_arr_ref, $percentage_thres);

my $major_allele_arr_ref=PositionProfiles::compute_mode_residue_for_profile($filtered_prof_arr_ref);
PositionProfiles::output_prof($output_prof_fname, $filtered_prof_arr_ref, $major_allele_arr_ref);

print STDERR "Done.\n";

