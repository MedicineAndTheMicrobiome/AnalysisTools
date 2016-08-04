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
use vars qw($opt_o $opt_r);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

my $clustalw_bin=`which clustalw2`;

getopts("o:r:");
my $usage = "usage: 
$0 

	-o <distFromRefence> -r <reference profile> <profile1> ... <profilen>

	Takes 1 to n profiles. 

	This script assumes that all your profiles are already pre-aligned, so 
	gaps already exist within and the LENGTH OF EACH PROFILE IS IDENTICAL.

	The output format is:

		<refName>\\t<prof1name>\\t<prof2name>\\t...\\t<profnname>\\n
		0.0000\\t<distance1,1>\\t<distance2,1>\\t...\\t<distancen,1>\\n
		0.0000\\t<distance1,2>\\t<distance2,2>\\t...\\t<distancen,2>\\n
		...
		0.0000\\t<distance1,length>\\t<distance2,length>\\t...\\t<distancen,length>\\n

	The first column is all zeros for historical reasons.
	It was a sanity check to ensure all profiles were compared against a reference.


";

if(!defined($opt_o) || !defined($opt_r)){
	die $usage;
}

my $output_file=$opt_o;
my $reference_prof=$opt_r;

my $prof_to_mfasta_bin=$FindBin::Bin;

###############################################################################

my @profile_names;

# Make sure everything is computed against the reference
print STDERR "Using $reference_prof as Reference.\n";
push @profile_names, $reference_prof;

# Include rest of the profiles
my $prof_name;
# Get profiles file names from command line
while($prof_name=shift){
	print STDERR "  Including $prof_name.\n";
	push @profile_names, $prof_name;
}

my $num_prof=$#profile_names+1;
print STDERR "Num profiles specified (including reference): $num_prof\n";

###############################################################################
# load profiles
my %loaded_major_seqs;
my %loaded_profiles; 

for(my $i=0; $i<$num_prof; $i++){
	print STDERR "Loading $profile_names[$i].\n";
	my ($seq_ref, $prof_ref)=PositionProfiles::load_prof($profile_names[$i]);
	$loaded_major_seqs{$i}=$seq_ref;
	$loaded_profiles{$i}=$prof_ref;
}

###############################################################################
# Compare against the reference(0)
my @results;
for(my $i=0; $i<=$#profile_names; $i++){
	print STDERR "Computing distance between $profile_names[$i] and Reference.\n";
	$results[$i]=PositionProfiles::compare_profiles_return_position_array("0", $i, \%loaded_profiles);
}

###############################################################################
# Write RMSD distance matrix out
print STDERR "Writing out results.\n";
open(OUT_FH, ">$output_file") || die "Could not open $output_file\n";

# Write header
my $header;
$header=join "\t", @profile_names;
print OUT_FH "$header\n";

# Write out positional distances
my $gapped_alignmen_length=$#{$results[0]}+1;
for(my $pos=0; $pos<$gapped_alignmen_length; $pos++){
	my @out=();
	for(my $i=0; $i<=$#profile_names; $i++){
		push @out, sprintf("%7.5f", ${$results[$i]}[$pos]);
	}
	my $outstr=join "\t", @out;
	print OUT_FH "$outstr\n";
}
close(OUT_FH);

###############################################################################

print STDERR "Done.\n";
