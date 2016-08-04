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
use vars qw($opt_o);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

my $clustalw2_bin=`which clustalw2`;
if($clustalw2_bin eq ""){
	die "\nCould not find clustalw2 program. Please specify in \$PATH environmental variable, or modifiy this script to help locate it.\n";
}

getopts("o:");
my $usage = "usage: 
$0 

	-o <output merged profile> <profile1> ... <profilen>

	Takes 1 to n profiles, and generates an output format with the help of
	clustalw.

	This program depends on an alignment program, such as clustalw.
	By default it is using: $clustalw2_bin

";

if(!defined($opt_o)){
	die $usage;
}

my $output_file=$opt_o;

my $prof_to_mfasta_bin=$FindBin::Bin;

###############################################################################

# Get profiles file names from command line
my @profile_names;
my $prof_name;
while($prof_name=shift){
	push @profile_names, $prof_name;
}
print STDERR "Num profiles specified: ", $#profile_names+1, "\n";

# Make temp directory to do all the work
my $tmp_dir="$output_file\.tmp";

if(-e $tmp_dir){
	print STDERR "Removing $tmp_dir to restart.\n";
	`rm -rf $tmp_dir`;
}

mkdir "$tmp_dir";

# Copy and count profiles to merge
my $num_prof=0;
foreach my $prof_name(@profile_names){
	`cp $prof_name $tmp_dir/$num_prof\.prof`;
	`chmod +w $tmp_dir/$num_prof\.prof`;
	$num_prof++;
}

# Convert profiles to multifasta for clustalw
print STDERR `$prof_to_mfasta_bin/Profiles_To_Multifasta.pl -o $tmp_dir/fasta $tmp_dir/*.prof`;

# Run clustalw
my $exec="$clustalw2_bin -infile=$tmp_dir/fasta -ktuple=4 -quicktree";
$exec=~s/\n//g;
print STDERR "Going to execute: $exec\n";
system($exec);

# Load clustalw output
my $aln_hash_ref=PositionProfiles::readClustalAlnFile("$tmp_dir/fasta.aln");

# load profiles
my %loaded_major_seqs;
my %loaded_profiles; 
for(my $i=0; $i<$num_prof; $i++){
	my ($seq_ref, $prof_ref)=PositionProfiles::load_prof($profile_names[$i]);
	$loaded_major_seqs{$i}=$seq_ref;
	$loaded_profiles{$i}=$prof_ref;
}

# Initialize blank profile
my $new_prof_length=$#{${$aln_hash_ref}{"1"}}+1;
print STDERR "New profile length: $new_prof_length\n";

my @new_prof_arr;

my ($gap_alphabet_ref, $nogap_alphabet_ref)=PositionProfiles::get_alphabet_profarr($loaded_profiles{0});
for(my $i=0; $i<$new_prof_length; $i++){
	foreach my $res(@{$gap_alphabet_ref}){
		${$new_prof_arr[$i]}{$res}=0;
	}
}

foreach my $seq_id (keys %{$aln_hash_ref}){
	print "Adding sequence $seq_id\n";
	PositionProfiles::add_profile(\@new_prof_arr, ${$aln_hash_ref}{$seq_id}, $loaded_profiles{$seq_id}, $loaded_major_seqs{$seq_id});
}

my $new_major_nuc_seq=PositionProfiles::compute_mode_residue_for_profile(\@new_prof_arr);

#print "$new_major_nuc_seq\n";

PositionProfiles::output_prof($output_file, \@new_prof_arr, $new_major_nuc_seq);

print STDERR "done.\n";

###############################################################################

