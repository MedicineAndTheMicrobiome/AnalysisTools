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
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Std;
use vars qw($opt_a $opt_g $opt_p);
use File::Basename;
use PositionProfiles;

getopts("a:g:d:p");
my $usage = "usage: 
$0 
	-a <NAST Alignment FASTA file>
	-g <Group Filename (1st col Sequence ID, 2nd col Group ID, tab separated)>
	[-p <use AA residues, default is nucleotide>]

	This script will read in a FASTA file that contains an aligned sequence
	using a NAST-based reference alignment.  A NAST alignment will have
	a bunch of gaps in it, becuase it was aligned to a reference, but each
	gap may not necessarily have the underlining residue that caused the gap
	to be opened, in this alignment file.  All sequences lengths will be the
	same (as a result of the NAST alignment).

	The group file contains a mapping of sequence id to group id.  There
	will be a profile build for each group id.

	Output will be a distance matrix of each group versus another group's
	RMSD.  
	
";

if(!defined($opt_a)||!defined($opt_g)){
	die $usage;
}

my $align_fasta=$opt_a;
my $group_file=$opt_g;

my $residue_type;
if(defined($opt_p)){
	$residue_type=PositionProfiles::AA_TYPE;
}else{
	$residue_type=PositionProfiles::NUC_TYPE;
}

print STDERR "Residue Type is: $residue_type\n";

###############################################################################

my %profiles;
my $ignored_sequences=0;

#------------------------------------------------------------------------------
# Load groups

my ($grp_hash_ref, $grp_arr_ref, $num_groups)=readGroupFile($group_file);

print STDERR "\n";
print STDERR "Num Groups: $num_groups\n";
print STDERR "Groups:\n";
for(my $i=0; $i<$num_groups; $i++){
	print("\t${$grp_arr_ref}[$i]\n");
}
print STDERR "\n";

#------------------------------------------------------------------------------
# Process alignment fasta

sub add_alignment_to_profile{
	my $profile_ref=shift;
	my $alignment_seq_ref=shift;

	my $seq=${$alignment_seq_ref};
	$seq=~s/\./-/g;

	my @residues=split //, $seq;
	
	# Remove from gaps
	for(my $i=0; $i<=$#residues; $i++){
		if($residues[$i] eq "-"){
			$residues[$i]=" ";
		}else{
			last;
		}
	}

	# Remove end gaps
	for(my $i=$#residues; $i>=0; $i--){
		if($residues[$i] eq "-"){
			$residues[$i]=" ";
		}else{
			last;
		}
	}

	# Accumulate counts into profile
	for(my $i=0; $i<=$#residues; $i++){
		if($residues[$i] ne " "){
			${${$profile_ref}[$i]}{$residues[$i]}++;
		}
	}
	
	return;
}

sub process_record{
	my $defline=shift;
	my $sequence=shift;

	$defline=~s/^>//;
	my ($seq_id)=split /\s+/, $defline;
	
	# Get Group ID
	my $grp_id=${$grp_hash_ref}{$seq_id};
	if(!defined($grp_id)){
		print STDERR "Ignored $seq_id.  Could not find group.\n";
		$ignored_sequences++;
		return;
	}

	#print "$defline\n";
	#print "$grp_id\n";
	#print "$sequence\n";

	if(!defined($profiles{$grp_id})){
		my @profile;
		my $seq_len=length($sequence);
		for(my $i=0; $i<$seq_len; $i++){
			$profile[$i]=PositionProfiles::init_pos_prof($residue_type);
		}
		$profiles{$grp_id}=\@profile;
		#print STDERR $profiles{$grp_id} . "\n";
		
	}

	add_alignment_to_profile($profiles{$grp_id}, \$sequence);
}

#------------------------------------------------------------------------------

open(FH, "<$align_fasta") || die "Could not open $align_fasta\n";

print STDERR "Processing FASTA file...\n";

my ($defline, $prev_defline, $sequence);
while(<FH>){
        chomp;

        if(/^>/){
                $defline=$_;
                if($sequence ne ""){
                        process_record($prev_defline, $sequence);
                        $sequence="";
                }
                $prev_defline=$defline;
        }else{
                $_=~s/\s+//g;
                $sequence.=$_;
        }
}
process_record($prev_defline, $sequence);

close(FH);
#------------------------------------------------------------------------------

foreach my $grp (sort @{$grp_arr_ref}){
	if(!defined($profiles{$grp})){
		print STDERR "No sequences for $grp.\n";
	}else{
		print STDERR "Writing profile for $grp.\n";
		my $profname="$grp\.prof";
		my $major_allele_seq=PositionProfiles::compute_mode_residue_for_profile($profiles{$grp});
		PositionProfiles::output_prof($profname, $profiles{$grp}, $major_allele_seq, $residue_type);
	}
}

#------------------------------------------------------------------------------


print STDERR "Ignored sequences: $ignored_sequences\n";
print STDERR "Done.\n";

###############################################################################
###############################################################################

sub readGroupFile{
	my $group_file=shift;
	print STDERR "Reading Group File: $group_file\n";
	open(GRP_FH, "<$group_file") || die "Could not open $group_file\n";

	my %group_hash;
	my %group_ids;

	while(<GRP_FH>){
		chomp;
		my ($seq_id, $grp_id)=split /\t/, $_;
		$group_hash{$seq_id}=$grp_id;
		$group_ids{$grp_id}=1;
	}

	my @group_id_arr=keys %group_ids;
	my $num_group_ids=$#group_id_arr+1;

	return(\%group_hash, \@group_id_arr, $num_group_ids);
}

###############################################################################
