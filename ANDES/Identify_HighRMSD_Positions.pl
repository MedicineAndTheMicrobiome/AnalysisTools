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
use vars qw($opt_r $opt_t $opt_o);

getopts("r:t:o:");
my $usage = "usage: 
$0 
	-r <RMSD output file>
	-t <threshold for filtering>
	-o <positions of interest file>

	The output will be:

	<position>\\<MembersMeetingCutff>\\t<RMSD1>\\t<RMSD2>\\t...
	
";

if(!defined($opt_r) || !defined($opt_r)){
	die $usage;
}

my $rmsd_in=$opt_r;
my $filt_out=$opt_o;
my $thresh=$opt_t;

open(IN_FH, "<$rmsd_in") || die "Could not open $rmsd_in\n";
open(OUT_FH, ">$filt_out") || die "Coult not open $filt_out\n";

my $first_line=<IN_FH>;
$first_line=~s/\n$//;
my @members=split /\t/, $first_line;
my @member_ids=split //, "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
my $reference=shift @members;

my $queries_str="";
for(my $i=0; $i<=$#members; $i++){
	$queries_str.="\t" . $member_ids[$i] . ":" . $members[$i];
}

print OUT_FH "#position\t$queries_str\n";

my $pos=0;
while(<IN_FH>){
	chomp;
	my @rmsds=split /\t/;
	shift @rmsds;

	my $idlist="";
	my $num_exc_thres=0;
	for(my $r=0; $r<=$#rmsds; $r++){
		if($rmsds[$r]>$thresh){
			$idlist.=$member_ids[$r];
			$num_exc_thres++;
		}
	}

	if($num_exc_thres>0){
		my $rmsd_str=join "\t", @rmsds;
		print OUT_FH "$pos\t$idlist\t$rmsd_str\n";
	}
	
	$pos++;
}


