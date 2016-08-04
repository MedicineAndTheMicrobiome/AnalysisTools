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
use vars qw($opt_i $opt_o $opt_t);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;
use Statistics::Descriptive;

getopts("i:o:t");
my $usage = "usage: 
$0 
	-i <Input Profile>
	-o <Output QVs>
	[-t type flag]

	Reads in the position profile and outputs the estimated QV values, assuming
	all the reads are supposed to be for the same sequence.

	OUTPUT:

		There will be two files generated:
			<Output QVs>.qv
			<Output QVs>.stat

		The .qv file will contain a QV for each position in your profile.  See below for 
		additional information if you set the -t flag.

		The .stat file will contain summary statistics of the qv's and a histogram.

	OTHER OPTIONS:

		If the -t option is specified, the second column in the .qv file will be one of
		the following:

			ins 	(Insertion error)
			del 	(Deletion error)
			sub	(Substitution error)
			con	(No error, complete Consensus at position)
			unk	(Unknown, possible error in profile)

		If the error is a substitution, then the 3rd column will contain the major
		allele, and the 4th column will contain the most common secondary/error allele:
		
		    24.0012 sub     A       G

		Means the QV is 24.0012, and it was a substitution error where an A was converted
		to a G.

		    24.0012 sub     G       Multi

		Means the G was converted to more than one alternate base.


		If the error is an indel:

		    26.9897 ins     1       A

		The majority of reads do not have the A allele, and the total number of A's up or downstream
		of that position is 1.
		
		    20.7831 del     2       G

		The majority of reads have a G.  The total up and downstream of that position are two G's.

";

if(
	!defined($opt_i) ||
	!defined($opt_o)
){
	die $usage;
}

my $input_prof_fname=$opt_i;
my $output_fname=$opt_o;
my $show_error_type=defined($opt_t)?1:0;

###############################################################################

my ($major_allele_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($input_prof_fname);

my $prof_length=$#{$prof_arr_ref}+1;
my @result_arr;

for(my $i=0; $i<$prof_length; $i++){
	$result_arr[$i]=PositionProfiles::estimate_qv(${$prof_arr_ref}[$i]);
}

# Compute some descriptive statistics
my $stat_obj=Statistics::Descriptive::Full->new();
$stat_obj->add_data(@result_arr);
my $mean=sprintf("%3.4f",$stat_obj->mean());
my $stdev=sprintf("%3.4f",$stat_obj->standard_deviation());
my $min=sprintf("%3.4f",$stat_obj->min());
my $max=sprintf("%3.4f",$stat_obj->max());
my %histogram=$stat_obj->frequency_distribution(30);
my $q1=sprintf("%3.4f",$stat_obj->percentile(25));
my $q2=sprintf("%3.4f",$stat_obj->percentile(50));
my $q3=sprintf("%3.4f",$stat_obj->percentile(75));
my $count=$stat_obj->count();

###############################################################################

my ($allele1_ref, $freq1_ref, $allele2_ref, $freq2_ref) = PositionProfiles::compute_top_allele_frequencies($prof_arr_ref);

my @error_type;
my @primary_allele;
my @secondary_allele;

my %error_type_table;
#error_type_table{primary}{secondary}="Status";

my $ERR_SUB="sub"; # Substitution Error
my $ERR_DEL="del"; # Deletion Error
my $ERR_INS="ins"; # Insertion Error
my $ERR_CON="con"; # Consensus
my $ERR_UNK="unk"; # Unknown Error

my $GAP="gap";       # All gaps
my $SINGLE="single"; # Single base
my $MULTI="multi";   # Multiple base
my $NULL="null";     # Null

$error_type_table{$SINGLE}{$SINGLE}=$ERR_SUB;
$error_type_table{$SINGLE}{$MULTI} =$ERR_SUB;
$error_type_table{$SINGLE}{$GAP}   =$ERR_DEL;
$error_type_table{$SINGLE}{$NULL}  =$ERR_CON;
$error_type_table{$MULTI }{$SINGLE}=$ERR_SUB;
$error_type_table{$MULTI }{$MULTI} =$ERR_SUB;
$error_type_table{$MULTI }{$GAP}   =$ERR_SUB;
$error_type_table{$MULTI }{$NULL}  =$ERR_SUB;
$error_type_table{$GAP   }{$SINGLE}=$ERR_INS;
$error_type_table{$GAP   }{$MULTI} =$ERR_INS;
$error_type_table{$GAP   }{$GAP}   =$ERR_UNK;
$error_type_table{$GAP   }{$NULL}  =$ERR_UNK;
$error_type_table{$NULL  }{$SINGLE}=$ERR_UNK;
$error_type_table{$NULL  }{$MULTI} =$ERR_UNK;
$error_type_table{$NULL  }{$GAP}   =$ERR_UNK;
$error_type_table{$NULL  }{$NULL}  =$ERR_UNK;

sub type{
	my $nuc=shift;
	my $freq=shift;

	my $type=undef;
	my $num_alleles=length($nuc);

	if(!($freq==0 || !defined($freq))){
		if($num_alleles==1){
			if($nuc eq "-"){
				$type=$GAP;
			}else{
				$type=$SINGLE;
			}
		}else{
			$type=$MULTI;
		}
	}else{
		$type=$NULL;
	}
	
	return($type);
}

my @error_type;
for(my $i=0; $i<$prof_length; $i++){
	
	my $pri_type=type(${$allele1_ref}[$i], ${$freq1_ref}[$i]);
	my $sec_type=type(${$allele2_ref}[$i], ${$freq2_ref}[$i]);

	my $error=$error_type_table{$pri_type}{$sec_type};
	#print "$pri_type / $sec_type\t\t$error\t\t";
	#PositionProfiles::print_pos_prof(${$prof_arr_ref}[$i]);

	my $error_info;
	my $allele_of_error;
	if($error eq $ERR_DEL || $error eq $ERR_INS){
			
		if($error eq $ERR_DEL){
			$allele_of_error=${$allele1_ref}[$i];
		}else{
			$allele_of_error=${$allele2_ref}[$i];
		}

		if(length($allele_of_error)==1){
			my $hp_length=PositionProfiles::is_indel_from_homopolymer($i, $allele_of_error, $prof_arr_ref, .5);
			$error_info.="$hp_length\t$allele_of_error";
		}else{
			$error_info.="\tMulti";
		}
	}elsif($error eq $ERR_SUB){
		my $major_allele=${$allele1_ref}[$i];
		my $allele_of_error=${$allele2_ref}[$i];
		if(length($major_allele)>1){
			$major_allele="Multi";	
		}
		if(length($allele_of_error)>1){
			$allele_of_error="Multi";
		}
		$error_info.="$major_allele\t$allele_of_error";
	}
	
	push @error_type, "$error\t$error_info";

}

###############################################################################

open(OUT_FH, ">$output_fname\.qv") || die "Could not open $output_fname\.qv for writing.\n";
my $error_out="";
for(my $i=0; $i<$prof_length; $i++){

	
	if($show_error_type){
		$error_out="\t$error_type[$i]";
	}

	print OUT_FH "$result_arr[$i]$error_out\n";
}
close(OUT_FH);

#-----------------------------------------------------------------------------

open(OUT_STAT_FH, ">$output_fname\.stat") || die "Could not open $output_fname\.stat for writing.\n";

print OUT_STAT_FH "\tmean\tstdev\tmin\tQ1\tQ2\tQ3\tmax\tcount\n";
print OUT_STAT_FH "Summary\t$mean\t$stdev\t$min\t$q1\t$q2\t$q3\t$max\t$count\n";
print OUT_STAT_FH "\n";
foreach my $key (sort {$a <=> $b} keys %histogram){
	print OUT_STAT_FH "$key\t$histogram{$key}\n";
}

close(OUT_STAT_FH);

###############################################################################

print STDERR "Done.\n";
