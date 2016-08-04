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
use BLOSUM_Utilities;

my $BU=new BLOSUM_Utilities("../BLOSUM/BLOSUM62");

$BU->print_matrix();
my $dist;
my $a;
my $b;

print "Scores:\n";

$a="ARND"; $b="ARND";
$dist=$BU->twoSequenceScore($a,$b);
print "$a/$b => $dist\n";

$a="ARNP"; $b="ARND";
$dist=$BU->twoSequenceScore($a,$b);
print "$a/$b => $dist\n";

$a="AR-D"; $b="ARND";
$dist=$BU->twoSequenceScore($a,$b);
print "$a/$b => $dist\n";

$a="HAPPYY"; $b="SADSAD";
$dist=$BU->twoSequenceScore($a,$b);
print "$a/$b => $dist\n";

$a="HAPPYY"; $b="SADSAD";
my @mask=(1,1,1,0,0,0);
$dist=$BU->twoSequenceScore($a,$b,\@mask);
print "$a/$b => $dist\n";


print "Dist:\n";

$a="ARND"; $b="ARND";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="ARNP"; $b="ARND";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="AR-D"; $b="ARND";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="HAPPYY"; $b="SADSAD";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="   CAT   "; $b="FATCATWOW";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="CAMPLETELY"; $b="DIFFERENT ";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="CAMPLETELY"; $b="DIFFERENT ";
my @mask=(1,1,1,1,1,0,0,0,0,0);
$dist=$BU->twoSequenceDist($a,$b, \@mask);
print "$a/$b => $dist\n";

$a="FF"; $b="YY";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="FFFF"; $b="YYYY";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";

$a="THIS-IS-HALLW  "; $b="THIS-IS-GWWDBYA";
my @mask=(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0);
$dist=$BU->twoSequenceDist($a,$b, \@mask);
print "$a/$b => $dist\n";

$a="THIS-IS-HALLW  "; $b="THIS-IS-GWWDBYA";
$dist=$BU->twoSequenceDist($a,$b);
print "$a/$b => $dist\n";


