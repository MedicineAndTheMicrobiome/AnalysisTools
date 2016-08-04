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
use vars qw($opt_o $opt_d);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

my $clustalw2_bin=`which clustalw2`;
if($clustalw2_bin eq ""){
        die "\nCould not find clustalw2 program. Please specify in \$PATH environmental variable, or modifiy this script to help locate it.\n";
}

getopts("o:d:");
my $usage = "usage: 
$0 

	-o <comparison matrix> <profile1> ... <profilen>
	[-d <level>]

	The -d option is useful if all your files have the same name, but are in
	different directories.  In this case, the name of the directory will be 
	used to label the output.  Level 1 is the directory immediate above the profile.
	Level 2 is the directory if you were to go two directories up.

	Takes 1 to n profiles, uses clustalw to align them, and then	
	computes the Root Mean Square Deviation (RMSD) between them.

	The input format is:

		<major variation>\\t<A count>\\t<T count>\\t<G count>\\t<C count>\\n

	The output format is:
		<Name1>\\t<Name2>\\t<Name3>\\t...\\t<Namen>\\n
		<score1,1>\\t<score2,1>\\t...\\t<score1,n>\\n
		<score1,2>\\t<score2,2>\\t...\\t<score2,n>\\n
		<score1,3>\\t<score2,3>\\t...\\t<score3,n>\\n
		...
		<score1,n>\\t<score2,n>\\t...\\t<scoren,n>\\n
	
	This program depends on an alignment program, such as clustalw2.
	By default it is using: $clustalw2_bin

";

if(!defined($opt_o)){
	die $usage;
}

my $use_dirname=undef;
if(defined($opt_d)){
	$use_dirname=$opt_d;	
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
	print STDERR `rm -f $tmp_dir/*.prof`;
	print STDERR `rm -f $tmp_dir/fasta`;
	print STDERR `rm -f $tmp_dir/fasta.aln`;
	print STDERR `rm -f $tmp_dir/fasta.dnd`;
}else{
	mkdir "$tmp_dir";
}

# Copy and count profiles to merge
my $num_prof=0;
my %prof_name_hash;
my %path_hash;
foreach my $prof_name(@profile_names){
	($prof_name_hash{$num_prof}, $path_hash{$num_prof})=fileparse($prof_name);

	if(defined($use_dirname)){
		$path_hash{$num_prof}=~s/\/$//;
		$path_hash{$num_prof}=~s/\/\//\//;
		my @path_components=split /\//, $path_hash{$num_prof};
		$prof_name_hash{$num_prof}=$path_components[($#path_components+1)-$use_dirname];
		print STDERR "Using directory as profile name: $prof_name_hash{$num_prof}\n";
	}

	$prof_name_hash{$num_prof}=~s/\.prof$//;

	my $cmd="cp $prof_name $tmp_dir/$num_prof\.prof";
	print STDERR "Going to execute: $cmd\n";
	print STDERR `$cmd`;

	#my $cmd="chmod +w $tmp_dir/$num_prof\.prof";
	#print STDERR "Going to execute: $cmd\n";
	#print STDERR `$cmd`;

	$num_prof++;
}

# Convert profiles to multifasta for clustalw
print STDERR "Converting profiles to major allele multifastas...\n";
if(`which Profiles_To_Multifasta.pl` eq ""){
	die "Could not find Profiles_To_Multifasta.pl.\n";
}
print STDERR `$prof_to_mfasta_bin/Profiles_To_Multifasta.pl -o $tmp_dir/fasta $tmp_dir/*.prof`;

# Run clustalw
print STDERR "Running MSA: $clustalw2_bin...\n";
#print STDERR `$clustalw2_bin -infile=$tmp_dir/fasta -quicktree -ktuple=4`;
my $exec="$clustalw2_bin -infile=$tmp_dir/fasta -quicktree -ktuple=4";
$exec=~s/\n//g;
print STDERR "Going to execute: $exec\n";
system($exec);

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
	foreach my $seq_id_b (keys %{$aln_hash_ref}){
		$results{$seq_id_a}{$seq_id_b}=PositionProfiles::compare_profiles($seq_id_a, $seq_id_b, \%gapped_profiles);
	}
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
foreach my $seq_id_a (sort keys %{$aln_hash_ref}){
	my @out=();
	foreach my $seq_id_b (sort keys %{$aln_hash_ref}){
		push @out, $results{$seq_id_a}{$seq_id_b};
	}
	my $outstr=join "\t", @out;
	print OUT_FH "$outstr\n";
}
close(OUT_FH);

#print STDERR `rm -rf $tmp_dir`;


