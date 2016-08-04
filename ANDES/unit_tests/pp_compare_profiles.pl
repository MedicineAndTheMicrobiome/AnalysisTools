#!/usr/bin/env perl

###############################################################################
#                                                                             # 
#       Copyright (c) 2005 J. Craig Venter Institute.                         #     
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

use FindBin ();
use lib "$FindBin::Bin/..";
use strict;
use PositionProfiles;


my ($major_nuc1_arr_ref, $prof1_arr_ref)=PositionProfiles::load_prof("short_nuc.prof");
my ($major_nuc2_arr_ref, $prof2_arr_ref)=PositionProfiles::load_prof("short_nuc2.prof");

my %prof_hash;
$prof_hash{"A"}=$prof1_arr_ref;
$prof_hash{"B"}=$prof2_arr_ref;

print "A\n";
PositionProfiles::print_prof($prof1_arr_ref);
print "B\n";
PositionProfiles::print_prof($prof2_arr_ref);

my $score=PositionProfiles::compare_profiles("A","A",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("A","B",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("B","A",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("B","B",\%prof_hash);
print "Score: $score\n";

my ($major_nuc1_arr_ref, $prof1_arr_ref)=PositionProfiles::load_prof("short_aa.prof");
my ($major_nuc2_arr_ref, $prof2_arr_ref)=PositionProfiles::load_prof("short_aa2.prof");

my %prof_hash;
$prof_hash{"A"}=$prof1_arr_ref;
$prof_hash{"B"}=$prof2_arr_ref;

print "A\n";
PositionProfiles::print_prof($prof1_arr_ref);
print "B\n";
PositionProfiles::print_prof($prof2_arr_ref);

my $score=PositionProfiles::compare_profiles("A","A",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("A","B",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("B","A",\%prof_hash);
print "Score: $score\n";
my $score=PositionProfiles::compare_profiles("B","B",\%prof_hash);
print "Score: $score\n";
