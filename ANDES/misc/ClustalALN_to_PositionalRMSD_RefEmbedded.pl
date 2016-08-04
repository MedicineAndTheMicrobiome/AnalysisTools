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
use vars qw($opt_r $opt_o);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("r:o:");
my $usage = "usage: 
$0 
	-r <reference prefix>
	-o <other prefix>
	<associated aln files...>

	Takes clustalw aln file, and then for each position, computes the rmsd
	between the reference set and the other set.

	Output is just a list of rmsd for every position.
	
";

if(!defined($opt_o) || !defined($opt_r)){
	die $usage;
}

my $ref_prefix=$opt_r;
my $other_prefix=$opt_o;

###############################################################################

my $aln_file;
my @fname_arr;

my $other_reads_profile_ref;
my $refer_reads_profile_ref;

my @unfiltered_profiles;

while($aln_file=shift){
	print STDERR "Loading $aln_file\n";
	$other_reads_profile_ref=PositionProfiles::read_ClustalALN_into_profile($aln_file, "$other_prefix");
	$refer_reads_profile_ref=PositionProfiles::read_ClustalALN_into_profile($aln_file, "$ref_prefix");
	push @unfiltered_profiles, $other_reads_profile_ref;
	push @fname_arr, $aln_file;
}

push @unfiltered_profiles, $refer_reads_profile_ref;
my $ref_prof_pos=$#unfiltered_profiles;

my @filtered_profiles=
	PositionProfiles::rebuild_filtered_all_gaps(@unfiltered_profiles);

my @rmsd_arr;
for(my $i=0; $i<=$ref_prof_pos; $i++){
	my $rmsd_position_arr_ref=PositionProfiles::compare_profiles_return_position_array_by_reference(
		    $filtered_profiles[$i], $filtered_profiles[$ref_prof_pos]);
	push @rmsd_arr, $rmsd_position_arr_ref;
}

print "reference";
for(my $j=0; $j<$ref_prof_pos; $j++){
	print "\t$fname_arr[$j]";
}
print "\n";

for(my $i=0; $i<=$#{$filtered_profiles[$ref_prof_pos]}; $i++){

	my $score = sprintf("%6.5f", ${$rmsd_arr[$ref_prof_pos]}[$i]);
        print "$score";

	for(my $j=0; $j<$ref_prof_pos; $j++){
		my $score = sprintf("%6.5f", ${$rmsd_arr[$j]}[$i]);
		print "\t$score";
	}
	print "\n";
}

print STDERR "Done!\n";

###############################################################################

