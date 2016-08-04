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
use FindBin;
use lib "$FindBin::Bin";
use PositionProfiles;

my $clustalw2_bin=`which clustalw2`;

getopts("o:");
my $usage = "usage: 
$0 

	-o <Output Conversion Table> <profile 1> <profile 2>

	Takes 2 profiles, uses clustalw to align them pairwise.
	Then generates a table per profile pair which contains the 
	conditional probability of residue conversions.

	This program uses an MSA program to align the profiles: $clustalw2_bin

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
mkdir "$tmp_dir";

# Copy and count profiles to merge
my $num_prof=0;
my %prof_name_hash;
foreach my $prof_name(@profile_names){
	($prof_name_hash{"$num_prof"})=fileparse($prof_name);
	$prof_name_hash{"$num_prof"}=~s/\.prof$//;
	print STDERR "[$num_prof] $prof_name_hash{$num_prof}\n";
	`cp $prof_name $tmp_dir/$num_prof\.prof`;
	`chmod +w $tmp_dir/$num_prof\.prof`;
	$num_prof++;
}

# Convert profiles to multifasta for clustalw
print STDERR `$prof_to_mfasta_bin/Profiles_To_Multifasta.pl -o $tmp_dir/fasta $tmp_dir/*.prof`;

# Run clustalw
my $cmd="$clustalw2_bin -infile=$tmp_dir/fasta";
$cmd=~s/\n//g;
print STDERR `$cmd`;

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

###############################################################################
# Compare all profiles versus all profiles

sub numerically{
	int($a) <=> int($b);
}

open(OUT_FH, ">$output_file") || die "Could not open $output_file\n";

my %trans_results;
my @id_arr=sort numerically keys %{$aln_hash_ref};

foreach my $seq_id_A (@id_arr){
	foreach my $seq_id_B (@id_arr){
		print STDERR "Computing Base Conversions between $seq_id_A and $seq_id_B.\n";
		$trans_results{$seq_id_A}{$seq_id_B}=PositionProfiles::compute_overall_transition_statistics($seq_id_A, $seq_id_B, \%gapped_profiles);
	}
}

my $first_seq_id=$id_arr[0];
my ($gap_res_arr_ref, $ungap_res_arr_ref)=PositionProfiles::get_alphabet_profarr($gapped_profiles{$first_seq_id});

#my @NUC_ARR=("A","T","G","C","-");
my $header=join "\t", @{$gap_res_arr_ref};

foreach my $seq_id_A (@id_arr){
        foreach my $seq_id_B (@id_arr){

		# Output header
		print OUT_FH $prof_name_hash{$seq_id_A} . "->" .  $prof_name_hash{$seq_id_B} . ": \n";
		print OUT_FH "\t$header\n";

		# Output the base change frequencies
		foreach my $src_nuc(@{$gap_res_arr_ref}){
			my @out;
			foreach my $dst_nuc(@{$gap_res_arr_ref}){
				my $nuc_trans_freq=${$trans_results{$seq_id_A}{$seq_id_B}}{$src_nuc}{$dst_nuc};
				push @out, sprintf("%5.2f%", $nuc_trans_freq*100.0);
			}
			my $outstr=join "\t", @out;
			print OUT_FH "$src_nuc\t$outstr\n";
		}

		print OUT_FH "\n";
	}
}

close(OUT_FH);

# Write transition & transversion probabilities, one per pair of nucleotides:
# Reminder:
#	Transition:	A <=> G  and  C <=> T
#	Transversion:	A <=> T,C  and  C <=> A,C
#	Indel:	A,T,G,C <=> gap

my %change_type;
$change_type{"A"}{"G"}="ts";
$change_type{"G"}{"A"}="ts";
$change_type{"C"}{"T"}="ts";
$change_type{"T"}{"C"}="ts";

###############################################################################

#print STDERR `rm -rf $tmp_dir`;


