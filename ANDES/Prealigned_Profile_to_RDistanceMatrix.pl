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
use vars qw($opt_p $opt_d $opt_n);
use File::Basename;
use PositionProfiles;

getopts("p:d:n");
my $usage = "usage: 
$0 
	-p <pre-aligned position profile list>
	-d <output distance matrix filename>
	[-n (don't trim off null positions)]

	Reads in a list of position profiles, then for each profile file name
	the profile is read in.

	If there are any null positions (all residues and gap are 0) among
	all profiles, they are trimmed off.  (This will increase the distances
	since, the distance is the average RMSD.)

	A distance matrix is then computed between all profiles and output.

";

if(!defined($opt_p)||!defined($opt_d)){
	die $usage;
}

my $profile_list=$opt_p;
my $dist_file=$opt_d;
my $trim_nulls;

if(defined($opt_n)){
	$trim_nulls=0;
}else{
	$trim_nulls=1;
}

###############################################################################

my $list_ref=readFileList($profile_list);
my @prof_arr=sort @{$list_ref};
my %profile_hash;

foreach my $profile_fn(@prof_arr){
	print STDERR "Loading $profile_fn\n";
	$profile_hash{$profile_fn}=PositionProfiles::load_prof($profile_fn);
	#PositionProfiles::print_prof($profile_hash{$profile_fn});
}
print STDERR "Completed loading all profiles.\n\n";

my $final_profiles_ref;
if($trim_nulls){

	#------------------------------------------------------------------------------
	# Determine how many ends to trim off

	my $prof_length=PositionProfiles::profile_arr_length($profile_hash{$prof_arr[0]});
	my $min_front_null=$prof_length;
	my $min_end_null=$prof_length;

	print STDERR "Determining how many positions on front and end of alignments may be trimmed.\n";
	foreach my $grp (@prof_arr){
		my ($grp_f_null, $grp_e_null)=PositionProfiles::size_null_ends($profile_hash{$grp});
		#print "$grp: Null Ends:  $grp_f_null / $grp_e_null\n";
		if($grp_f_null<$min_front_null){
			$min_front_null=$grp_f_null;
		}
		if($grp_e_null<$min_end_null){
			$min_end_null=$grp_e_null;
		}
	}

	print STDERR "Null ends to trim: From Beginning: $min_front_null / From End: $min_end_null\n\n";

	#------------------------------------------------------------------------------
	# Trim off ends

	my %trimmed_profiles;
	my $start=$min_front_null;
	my $stop=$prof_length-$min_end_null-1;

	print STDERR "Keeping from $start to $stop (inclusive).\n";

	foreach my $grp (@prof_arr){
		print STDERR "Trimming Null Ends: $grp\n";
		$trimmed_profiles{$grp}=PositionProfiles::normalize_profile(
			PositionProfiles::trim($profile_hash{$grp}, $start, $stop));
	}

	print STDERR "Trimming completed.\n\n";

	$final_profiles_ref=\%trimmed_profiles;
}else{
	$final_profiles_ref=\%profile_hash;
}

#------------------------------------------------------------------------------
# Compute distance matrix now:

my %dist_mat;

foreach my $grp1 (@prof_arr){
	foreach my $grp2 (@prof_arr){
		print STDERR "Comparing $grp1 vs $grp2.\n";
		$dist_mat{$grp1}{$grp2}=PositionProfiles::compare_profiles(
			$grp1, $grp2, $final_profiles_ref, 1);
	}
}

print STDERR "Comparisons complete.\n\n";

#------------------------------------------------------------------------------
# Output distance matrix

print STDERR "Outputing Distance Matrix.\n\n";
printMatrix(\%dist_mat, $dist_file);

print STDERR "\n\nDone.\n";

###############################################################################
###############################################################################

sub printMatrix{
	my $matrix_hash_ref=shift;
	my $outputfilename=shift;
	my $min_dist=shift;

	open(MAT_OUT, ">$outputfilename") || die "Could not open $outputfilename for writing.\n";
	my @matrix_members=sort keys %{$matrix_hash_ref};

	# print header
	my $headers=join " ", @matrix_members;
	print MAT_OUT " $headers\n";

	# print matrix info
	foreach my $a_id (@matrix_members){
		print MAT_OUT "$a_id";
		foreach my $b_id (@matrix_members){
			my $out_val=sprintf "%3.6f", ${$matrix_hash_ref}{$a_id}{$b_id};
			if($out_val <= 0 && $a_id ne $b_id){
				$out_val = sprintf "%3.6f", $min_dist;
			}
			print MAT_OUT " $out_val";
		}
		print MAT_OUT "\n";
	}

	close(MAT_OUT);
}

######################################################################################

sub readFileList{
	my $list_fn=shift;
	print STDERR "Reading List: $list_fn\n";
	open(LIST_FH, "<$list_fn") || die "Could not open $list_fn\n";

	my @list;

	while(<LIST_FH>){
		chomp;
		push @list, $_;
	}

	return(\@list);
}

###############################################################################
