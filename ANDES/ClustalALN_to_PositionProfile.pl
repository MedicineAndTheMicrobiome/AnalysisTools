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
use vars qw($opt_a $opt_m);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;


my $THRESHOLD=100;
my $MAJOR_FASTA="major.fasta";
my $PROFILE="prof";

getopts("a:m");
my $usage = "usage: 
$0 
	-a <.ALN file>
	[-m (aMino acid flag)]

	Takes aln file from clustalw and generates a profile file and major variation fasta file.

	By default, nucleotides are assumed.
	Use the -m, aMino acid flag, to specify using the amino acid alphabet.

	The output format of the profile file is:

	    <major variation>\\t<residue1 count>\\t<residue2 count\\t...<residueN count>\\t<gap count>\\n

	The residues are in alphabetical order. For example, for nucleotides, A, C, G, T.

	The filenames are:

	<aln file>.<type>.$MAJOR_FASTA
	<aln file>.<type>.$PROFILE

	Where <type>=aa for amino acids, or <type>=nuc for nucleotides.

	
";

if(!defined($opt_a)){
	die $usage;
}

my $aln_file=$opt_a;

my $residue_type=PositionProfiles::NUC_TYPE;
if(defined($opt_m)){
	$residue_type=PositionProfiles::AA_TYPE;
}

print STDERR "Residue Type: $residue_type\n";

my $fasta_output_file=$aln_file;
$fasta_output_file=~s/aln$/$residue_type\.$MAJOR_FASTA/;

my $profile_output_file=$aln_file;
$profile_output_file=~s/aln$/$residue_type\.$PROFILE/;

# Get sequence name for defline based on alignment file name
my ($seqname, $path)=fileparse($aln_file);
$seqname=~s/\.aln//;

###############################################################################

my $seq_hash_ref=PositionProfiles::readClustalAlnFile($aln_file);

###############################################################################

print STDERR "Summing Up Profiles...\n";
my $profile_ref=PositionProfiles::sequence_hash_into_profile($seq_hash_ref, $residue_type);

print STDERR "Computing Mode Residue...\n";
my $major_sequence=PositionProfiles::compute_mode_residue_for_profile($profile_ref);
	
###############################################################################

print STDERR "Writing Major Variations: $fasta_output_file\n";
PositionProfiles::output_fasta($fasta_output_file, $seqname, $major_sequence, 60);

print STDERR "Writing Profile: $profile_output_file\n";
PositionProfiles::output_prof($profile_output_file, $profile_ref, $major_sequence, $residue_type);

print STDERR "Done!\n";

###############################################################################
