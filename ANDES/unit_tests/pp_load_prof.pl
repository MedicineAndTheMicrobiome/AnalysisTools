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


my ($major_nuc_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof("short_nuc.prof");

my $str=join "", @{$major_nuc_arr_ref};
print "$str\n";

PositionProfiles::print_prof($prof_arr_ref);

my ($major_nuc_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof("short_aa.prof");

my $str=join "", @{$major_nuc_arr_ref};
print "$str\n";

PositionProfiles::print_prof($prof_arr_ref);
