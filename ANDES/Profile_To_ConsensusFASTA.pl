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
use vars qw($opt_c $opt_p $opt_g $opt_f);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("c:p:g:f:");
my $usage = "usage: 
$0 
	-c <output consensus fasta file> 
	-p <input profile>

	Fill gaps with N, instead of underlying residue.  This is essentially a placehold for gaps.
	E.g., .05 will change a position to an N if greater than 5% of the position is gap>]
	[-g <Threshold for substituting gaps with N>]

	Remove position if dominated with gaps. (ie. remove low frequency insertions)
	E.g., .95, will remove a position if there is only a 5% insertion rate>]
	[-f <Threshold for removing gaps.>] 

	The default is to take the underlying residues whatever the proportion of nucleotides to gap.
	So for example, if 95% gap, 2% G, and 3% T, the result for that position is K.

	Reads in a profile and generates a consensus FASTA file from it.
	The consensus fasta will have IUPAC ambiguity codes in it to
	represent degenerate bases.  Make sure you apply some filtering
	to the profiles ahead of time, or else you may get too many ambiguity
	codes.  Even if there is a very small representation of an alternate
	base, all non-zero frequency bases will be represented.
	
";

if(
	!defined($opt_c) ||
	!defined($opt_p)
){
	die $usage;
}

my $in_prof=$opt_p;
my $out_fasta=$opt_c;

my $gap_to_N_threshold=-1;
my $gap_removal_threshold=-1;

if(defined($opt_g)){
	$gap_to_N_threshold=$opt_g;
}

if(defined($opt_f)){
	$gap_removal_threshold=$opt_f;
}

if(defined($opt_g) && defined($opt_f)){
	die "Error:  You can not specify the gap thresholds for N substitution and remove simultaneously.\n";
}

print STDERR "Gap to N threshold: $gap_to_N_threshold\n";
print STDERR "Gap Filter: $gap_removal_threshold\n";

###############################################################################

my ($major_nuc_ref, $prof_ref)=PositionProfiles::load_prof($in_prof);

my $prof_length=$#{$prof_ref} +1;

print STDERR "Profile length = $prof_length\n";

my $consensus_seq="";

for(my $i=0; $i<$prof_length; $i++){

	my $normalized_prof_ref=PositionProfiles::normalize(${$prof_ref}[$i]);

	# If profile pos is null, ie. no residues or gaps defined then skip it.
	if(!PositionProfiles::is_null_prof($normalized_prof_ref)){

		my $perc_gap=${$normalized_prof_ref}{"-"};
	
		if($perc_gap!=1){

			my $consensus_nuc;
			if($gap_to_N_threshold!=-1 && $perc_gap>=$gap_to_N_threshold){
				$consensus_nuc="N";
			}elsif($gap_removal_threshold!=-1 && $perc_gap>=$gap_removal_threshold){
				$consensus_nuc="";
			}else{
				$consensus_nuc=PositionProfiles::prof_to_ambig(${$prof_ref}[$i]);
			}

			$consensus_seq.=$consensus_nuc;
		}
	}	
}

print STDERR "\n\nConsensus Length = " . length($consensus_seq) . "\n";

my ($defline)=fileparse($in_prof);
$defline=~s/\.prof$//;

output_fasta($out_fasta, $defline, $consensus_seq);

###############################################################################

sub output_fasta{
	my $fasta_fn=shift;
	my $defline=shift;
	my $sequence=shift;

	open(OUT_FH, ">$fasta_fn") || die "Could not open output fasta: $fasta_fn\n";

	print OUT_FH ">$defline\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print OUT_FH substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

	close(FH);

}
