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
chomp $clustalw_bin;

getopts("o:r:");
my $usage = "usage: 
$0 

	-o <distFromRefence> -r <reference profile> <profile1> ... <profilen>

	Takes 1 to n profiles, uses clustalw2 to align them, and then	
	computes the RMSD relative to the references.

	The input format is:

		<major variation>\\t<A count>\\t<T count>\\t<G count>\\t<C count>\\n

	The output format is:

		<refName>\\t<prof1name>\\t<prof2name>\\t...\\t<profnname>\\n
		0.0000\\t<distance1,1>\\t<distance2,1>\\t...\\t<distancen,1>\\n
		0.0000\\t<distance1,2>\\t<distance2,2>\\t...\\t<distancen,2>\\n
		...
		0.0000\\t<distance1,length>\\t<distance2,length>\\t...\\t<distancen,length>\\n

	The first column is all zeros for historical reasons.
	It was a sanity check to ensure all profiles were compared against a reference.
	
	Using clustalw2 from $clustalw_bin
	

";

if(!defined($opt_o) || !defined($opt_r)){
	die $usage;
}

my $output_file=$opt_o;
my $reference_prof=$opt_r;

my $prof_to_mfasta_bin=$FindBin::Bin;

###############################################################################

# Get profiles file names from command line
my @profile_names;

# Make sure everything is computed against the reference
push @profile_names, $reference_prof;

my $prof_name;
while($prof_name=shift){
	push @profile_names, $prof_name;
}
print STDERR "Num profiles specified (including reference): ", $#profile_names+1, "\n";

# Make temp directory to do all the work
my $tmp_dir="$output_file\.tmp";
mkdir "$tmp_dir";

# Copy and count profiles to merge
my $num_prof=0;
my %prof_name_hash;
foreach my $prof_name(@profile_names){
	($prof_name_hash{$num_prof})=fileparse($prof_name);
	$prof_name_hash{$num_prof}=~s/\.prof$//;
	print STDERR "[$num_prof] $prof_name_hash{$num_prof}\n";
	`cp $prof_name $tmp_dir/$num_prof\.prof`;
	`chmod +w $tmp_dir/$num_prof\.prof`;
	$num_prof++;
}

# Convert profiles to multifasta for clustalw
print STDERR "Converting profiles to multifasta\n";
print STDERR `$prof_to_mfasta_bin/Profiles_To_Multifasta.pl -o $tmp_dir/fasta $tmp_dir/*.prof`;

# Run clustalw
print STDERR "Running $clustalw_bin to generate alignments.\n";
print STDERR `$clustalw_bin -infile=$tmp_dir/fasta -quicktree`;
print STDERR "Finished with major allele alignments.\n";

# Load clustalw output
my $aln_hash_ref=PositionProfiles::readClustalAlnFile("$tmp_dir/fasta.aln");

#foreach my $id(keys %{$aln_hash_ref}){
#	my $str=join "", @{${$aln_hash_ref}{$id}};
#	#print "$id: $str\n";
#}

# load profiles
my %loaded_major_seqs;
my %loaded_profiles; 
for(my $i=0; $i<$num_prof; $i++){
	my ($seq_ref, $prof_ref)=PositionProfiles::load_prof($profile_names[$i]);
	$loaded_major_seqs{$i}=$seq_ref;
	$loaded_profiles{$i}=$prof_ref;
}

# Add gaps to profiles that need gaps
my %gapped_profiles;
foreach my $seq_id (keys %{$aln_hash_ref}){
	print STDERR "Adding gaps to sequence: $seq_id\n";
	$gapped_profiles{$seq_id}=
		PositionProfiles::add_gaps_to_profile(${$aln_hash_ref}{$seq_id}, $loaded_profiles{$seq_id});
}

# Compare all the gapped profiles against each other.  The results should be
#	symmetric, but we will double compute them for a sanity check.
my %results;
foreach my $seq_id_a (keys %{$aln_hash_ref}){
	$results{$seq_id_a}=PositionProfiles::compare_profiles_return_position_array("0", $seq_id_a, \%gapped_profiles);
}

# Write RMSD distance matrix out
open(OUT_FH, ">$output_file") || die "Could not open $output_file\n";

my @headers;
my $header;
foreach my $id (sort keys %{$aln_hash_ref}){
	push @headers, $prof_name_hash{$id};
}
$header=join "\t", @headers;

print OUT_FH "$header\n";
my $gapped_alignmen_length=$#{$results{"0"}}+1;


print "Gapped alignment length: $gapped_alignmen_length\n";
my @profiles=sort keys %{$aln_hash_ref};
my $profiles_str=join "/", @profiles;
print "Profiles: $profiles_str\n";

for(my $pos=0; $pos<$gapped_alignmen_length; $pos++){
	my @out=();
	foreach my $seq_id_a (@profiles){
		push @out, sprintf("%7.5f", ${$results{$seq_id_a}}[$pos]);
	}
	my $outstr=join "\t", @out;
	print OUT_FH "$outstr\n";
}
close(OUT_FH);

#print STDERR `rm -rf $tmp_dir`;


