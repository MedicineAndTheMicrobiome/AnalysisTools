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
use vars qw($opt_i);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("i:");
my $usage = "usage: 
$0 
	-i <Input Profile>

	Reads in the position profile and generates table like the following:

		[0]     [1]     [2]     [3]     [4]     [total]
	[A]     13      8       5       3       1       = 30
	[T]     7       9       5               1       = 22
	[G]     18      20      4       3               = 45
	[C]     5       7       5                       = 17

	Identifies deep gaps, and then the length and nucleotide of the homopolymer
	that probably caused it.  Then bins them up. 
	
";

if(
	!defined($opt_i)
){
	die $usage;
}

my $input_prof_fname=$opt_i;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

my $prof_length=$#{$prof_arr_ref}+1;

# Detect gaps
my $threshold=.95;
my $nucleotide="-";
my $normalized=1;
my $gap_pos_array_ref=PositionProfiles::find_residue_pos_above_threshold($prof_arr_ref, $threshold, $nucleotide);

# Collapse continuous gaps together
my $prev_pos;
my @collapsed_arr;
for(my $i=0; $i<($prof_length-1); $i++ ){
	if(${$gap_pos_array_ref}[$i]!=${$gap_pos_array_ref}[$i+1]){
		push @collapsed_arr, ${$gap_pos_array_ref}[$i];
	}
}
$gap_pos_array_ref=\@collapsed_arr;

# Step through the identified and collapse gaps
my $search_threshold=.05;
my %gap_attribute_hash;
my $max_stretch=0;
foreach my $pos(@{$gap_pos_array_ref}){

	# Get stats on gap position
	my $nuc_that_caused_gap = ${$major_allele_arr_ref}[$pos];
	my $freq=PositionProfiles::freq_at_pos($prof_arr_ref, $pos, $nuc_that_caused_gap, $normalized);
	my $nuc_freq=PositionProfiles::freq_at_pos($prof_arr_ref, $pos, $nuc_that_caused_gap, 0);

	#print "$pos: $nuc_that_caused_gap ($freq)\n";
	#print "Count at pos: $nuc_freq\n";

	# Count up and downstream nucleoties that are the same as the nuc that caused the gap	
	my $track_length=0;

	# Count reverse threshold
	my $spos=$pos-1;
	my $nuc_prob;
	do{
		$nuc_prob=PositionProfiles::freq_at_pos($prof_arr_ref, $spos, $nuc_that_caused_gap, 1);
		#print "\t$nuc_prob $nuc_that_caused_gap $spos\n";
		$track_length++;
		$spos--;
	}while($spos>=0 && $nuc_prob>=$search_threshold);

	#print "\t--\n";

	# Count forward threshold
	my $spos=$pos+1;
	my $nuc_prob;
	do{
		$nuc_prob=PositionProfiles::freq_at_pos($prof_arr_ref, $spos, $nuc_that_caused_gap, 1);
		#print "\t$nuc_prob $nuc_that_caused_gap $spos\n";
		$track_length++;
		$spos++;
	}while($spos<$prof_length && $nuc_prob>=$search_threshold);

	# Correct for adding too many in do/while loop
	$track_length-=2; 
	#print "    TrLen=$track_length\n";
	if($track_length>$max_stretch){
		$max_stretch=$track_length;
	}

	# Track frequency of track lengths for each nuc
	$gap_attribute_hash{$nuc_that_caused_gap}{$track_length}+=$nuc_freq;	

}

################################################################################

# For formatting table, if stretches didn't reach 10, pad the table out.
if($max_stretch<10){
	$max_stretch=10;
}

# Output headers
print "\n";
print "\nRaw counts:\n";
for (my $i=0; $i<=$max_stretch; $i++){
	print "\t[$i]";
}
print "\t[total]\n";
# Output frequencies
my %length_hash;
foreach my $nuc("A","T","G","C"){
	print "[$nuc]\t";
	my $nuc_total=0;
	for(my $i=0; $i<=$max_stretch; $i++){
		print "$gap_attribute_hash{$nuc}{$i}\t";
		$nuc_total+=$gap_attribute_hash{$nuc}{$i};
		$length_hash{$i}+=$gap_attribute_hash{$nuc}{$i};
		$length_hash{"total"}+=$gap_attribute_hash{$nuc}{$i};
	}
	print "= $nuc_total";
	print "\n";
}

print "\n[Total]";
for(my $i=0; $i<=$max_stretch; $i++){
	print "\t$length_hash{$i}";
}
print "\t= " . $length_hash{"total"};
print "\n\n";


# Compute nuc composition
my $nuc_composition=PositionProfiles::compute_residue_composition($prof_arr_ref, 0, $prof_length, .5);
print "\nNucleotide composition:\n";
foreach my $nuc("A","T","G","C"){
	print "\t$nuc: ${$nuc_composition}{$nuc}\n";
}

# Compute bases sequenced
my %nuc_probs;
my ($total_sequenced, $bases_sequenced_hash_ref)=PositionProfiles::compute_total_residues_assayed($prof_arr_ref, 0, $prof_length);
print "\nTotal bases sequenced: $total_sequenced\n";
foreach my $nuc("A","T","G","C"){
	my $perc= ${$bases_sequenced_hash_ref}{$nuc}/$total_sequenced;
	$nuc_probs{$nuc}=$perc;
	print "\t$nuc: ${$bases_sequenced_hash_ref}{$nuc} [$perc]\n";
}

# Normalized 
my $FACTOR=10000;
print "\nNormalized gaps per $FACTOR:\n";
for (my $i=0; $i<=$max_stretch; $i++){
	print "\t[$i]";
}
print "\t[total]\n";
# Output frequencies
foreach my $nuc("A","T","G","C"){
	print "NG [$nuc]\t";
	my $nuc_total=0;
	for(my $i=0; $i<=$max_stretch; $i++){
		my $norm=sprintf("%3.3f", $FACTOR*$gap_attribute_hash{$nuc}{$i}/${$bases_sequenced_hash_ref}{$nuc});
		print "$norm\t";
		$nuc_total+=$norm;
	}
	print "= $nuc_total";
	print "\n";
}

print STDERR "\nDone.\n";
