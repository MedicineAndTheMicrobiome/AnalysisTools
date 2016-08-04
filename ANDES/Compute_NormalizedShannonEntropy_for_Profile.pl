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
use vars qw($opt_i $opt_o $opt_g);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;
use Statistics::Descriptive;

getopts("i:o:g:");
my $usage = "usage: 
$0 
	-i <Input Profile>
	-o <Output Normalized Shannon Entropy Values root filename>
	[-g <Ignore positions with gaps exceeding threshold, eg .75>]

	Reads in the position profile and outputs the Shannon Entropy calculations
	for each position.

	The -g option forces positions which contain greater than specified threshold
	of gaps to be ignored.  So, if you set the threshold to .75, than if >75%
	of the position are gaps, ignore that position.
	
";

if(
	!defined($opt_i) ||
	!defined($opt_o)
){
	die $usage;
}

my $input_prof_fname=$opt_i;
my $output_entropy_fname=$opt_o;
my $gap_threshold=$opt_g;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

# Perform gap filtering
if(defined($gap_threshold)){
	if($gap_threshold >=0 && $gap_threshold<=1){
		$prof_arr_ref=PositionProfiles::filter_gaps($prof_arr_ref, $gap_threshold);
	}else{
		die "Gap threshold is out of range: $gap_threshold is not between 0 and 1\n";
	}
}

# Compute normalized shannon entropy
my $prof_length=$#{$prof_arr_ref}+1;
my @se_arr;

for(my $i=0; $i<$prof_length; $i++){
	$se_arr[$i]=PositionProfiles::normalized_shannon_entropy(${$prof_arr_ref}[$i]);
}

# Compute some descriptive statistics
my $stat_obj=Statistics::Descriptive::Full->new();
$stat_obj->add_data(@se_arr);
my $mean=sprintf("%3.4f",$stat_obj->mean());
my $stdev=sprintf("%3.4f",$stat_obj->standard_deviation());
my $min=sprintf("%3.4f",$stat_obj->min());
my $max=sprintf("%3.4f",$stat_obj->max());
my %histogram=$stat_obj->frequency_distribution(30);
my $q1=sprintf("%3.4f",$stat_obj->percentile(25));
my $q2=sprintf("%3.4f",$stat_obj->percentile(50));
my $q3=sprintf("%3.4f",$stat_obj->percentile(75));
my $count=$stat_obj->count();

###############################################################################

open(OUT_FH, ">$output_entropy_fname\.ent") || die "Could not open $output_entropy_fname\.ent for writing.\n";
for(my $i=0; $i<$prof_length; $i++){
	print OUT_FH "$i\t$se_arr[$i]\n";
}
close(OUT_FH);

#-----------------------------------------------------------------------------

open(OUT_STAT_FH, ">$output_entropy_fname\.stat") || die "Could not open $output_entropy_fname\.stat for writing.\n";

print OUT_STAT_FH "\tmean\tstdev\tmin\tQ1\tQ2\tQ3\tmax\tcount\n";
print OUT_STAT_FH "Summary\t$mean\t$stdev\t$min\t$q1\t$q2\t$q3\t$max\t$count\n";
print OUT_STAT_FH "\n";
foreach my $key (sort {$a <=> $b} keys %histogram){
	print OUT_STAT_FH "$key\t$histogram{$key}\n";
}

close(OUT_STAT_FH);

###############################################################################

print STDERR "Done.\n";
