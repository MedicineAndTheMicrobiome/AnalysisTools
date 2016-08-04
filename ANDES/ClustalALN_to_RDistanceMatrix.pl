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
use vars qw($opt_a $opt_s $opt_r $opt_m $opt_o $opt_M $opt_P $opt_g $opt_b $opt_c);
use File::Basename;
use BLOSUM_Utilities;

getopts("a:sr:m:o:M:Pg:b:c");
my $usage = "usage: 
$0 
	-a <.ALN file>
	[-P (Amino Acid/Protein Flag)]
	[-m <minimum non self-self distance]
	[-o <output file name>]
	
	Mask related options:
	[-M <Mask file name>:<Mask sequence id in .ALN file>]
	[-g <gap mask fill rule, \"zero\" or \"average\", default=average>]
	[-c (flag to not output Mask sequence id in distance matrix)]
	
	Nucleotide only options:
	[-s ignore stutter flag, don't treat gaps next to nucleotides that repeat as a difference.]
	[-r <min representation of variation to not be considered sequencing error, eg 5%>]

	Setting a different similarity (eg. BLOSUM) matrix
	[-b <peptide substitution matrix file name>]

	Takes aln file from clustalw and generates a distance matrix that can 
	be used in R.

	The distance is defined as (1-(percent identical)).

	N's match anything.

";

if(!defined($opt_a)){
	die $usage;
}

my $ignore_stutter=$opt_s;
my $aln_file=$opt_a;
my $output_file=$aln_file;
my $min_dist=0;
my $min_rep=defined($opt_r)?$opt_r:0;
my $dont_output_masksequenceid=defined($opt_c)?1:0;
my $is_AA_alignment=0;
my $gap_mask_fill_method="average";

my $peptide_similarity_matrix="$FindBin::Bin/BLOSUM/BLOSUM62_J";
if(defined($opt_b)){
	$peptide_similarity_matrix=$opt_b;
	print STDERR "Looking for $peptide_similarity_matrix\n";
	if(!(-e $peptide_similarity_matrix)){
		$peptide_similarity_matrix="$FindBin::Bin/BLOSUM/$peptide_similarity_matrix";
		print STDERR "Looking for $peptide_similarity_matrix\n";
		if(!(-e $peptide_similarity_matrix)){
			die "Could not find $opt_b\n";
		}
	}
}

my $BU=new BLOSUM_Utilities($peptide_similarity_matrix);

if(defined($opt_m)){
	$min_dist=$opt_m;
}

if(defined($opt_P)){
	$is_AA_alignment=1;
}

if(defined($opt_o)){
	$output_file=$opt_o;
}else{
	$output_file=~s/aln$/r_distmat/;
}

if(defined($opt_g)){
	$gap_mask_fill_method=$opt_g;
}

#-------------------------------------------------------------------------------
# Mask related setup
my $mask_supplied=0;
my $mask_fname;
my $mask_seqname;
if(defined($opt_M)){
	($mask_fname, $mask_seqname)=split /:/, $opt_M;	
	if($mask_seqname eq "" || !defined($mask_seqname)){
		$mask_seqname="mask_reference";
		print STDERR "Undefined mask sequence, using $mask_seqname.\n";
	}
	print STDERR "Mask Filename: $mask_fname\n";
	print STDERR "Mask Sequence Name: $mask_seqname\n";
	$mask_supplied=1;
}else{
	print STDERR "Mask not defined.\n";
}
#-------------------------------------------------------------------------------

print STDERR "Output is going to: $output_file\n";

my ($seqname, $path)=fileparse($aln_file);
$seqname=~s/\.aln//;

###############################################################################

# Load mask values if they were supplied
my $mask_val_arr_ref;
if($mask_supplied){
	$mask_val_arr_ref=load_mask($mask_fname);
	#print_mask_values($mask_val_arr_ref);
}else{

}

###############################################################################

print STDERR "Reading ALN file...\n";

open(ALN_FH, "<$aln_file") || die "Could not open $aln_file\n";

my %seq_hash;
my $clustal_file=0;
while(<ALN_FH>){
	chomp;
	if($_=~/CLUSTAL/){
		$clustal_file=1;
		next;	
	}elsif($_=~/\*/){
		next;
	}else{
		my ($id, $sequence)=split /\s+/, $_;
		if($id ne ""){
			$seq_hash{$id}.=uc($sequence);
		}
	}
}
close(ALN_FH);

if($clustal_file!=1){
	print STDERR "Are you sure this is clustal file?\n";
}

print STDERR "Done...\n";

###############################################################################


my %distance_hash;

# Pre change all sequence strings into arrays after remove leading/trailing gaps
my $seq_len=-1;
my $amb_residue;

# Select the ambiguous base based on the residue type
if($is_AA_alignment){
	$amb_residue=" ";
}else{
	$amb_residue="N";
}

# For each sequence, clean it up.
foreach my $id (keys %seq_hash){
	# Remove leading gaps
	if($seq_hash{$id}=~s/^(-+)//){
		$seq_hash{$id}=($amb_residue x length($1)) . $seq_hash{$id};
	}

	# Remove trailing gaps
	if($seq_hash{$id}=~s/(-+)$//){
		$seq_hash{$id}.=($amb_residue x length($1));
	}

	# Convert sequence strings to arrs
	my @arr=split //, $seq_hash{$id}; #/
	$seq_hash{$id}=\@arr;

	# Make sure all lengths are the same
	if($seq_len==-1){
		$seq_len=$#arr+1;
	}else{
		if($seq_len!=$#arr+1){
			die "Length of all sequences do not appear to be the same.\n";
		}
	}
}

#-----------------------------------------------------------------------------------

if(!$is_AA_alignment){

	# Compute array of acceptable nucleotides
	my @rep_arr; # for example: ${$rep_arr[$n]}{"A"}=1
	for(my $i=0; $i<$seq_len; $i++){
		my %comp_hash;
		my $totals=0;

		# Count up the number of representatives for each variation
		foreach my $seq_id(keys %seq_hash){
			my $nuc=${$seq_hash{$seq_id}}[$i];
			if(!defined($comp_hash{$nuc})){
				$comp_hash{$nuc}=0;
			}
			$comp_hash{$nuc}++;
			$totals++;
		}
		
		# See if the number of representatives is greater than the threshold, if so 
		#  then assume it's not a sequencing error.
		foreach my $nuc(keys %comp_hash){
			if(100.0*$comp_hash{$nuc}/$totals > $min_rep){
				${$rep_arr[$i]}{$nuc}=1;
			}
		}

	}

	# Replace nucletides that may be sequencing error with spaces, to differentiate from gaps/-'s
	foreach my $id (keys %seq_hash){
		for(my $i=0; $i<$seq_len; $i++){
			my $nuc=${$seq_hash{$id}}[$i];
			if($nuc=~/[CATG]/){
				if(!defined(${$rep_arr[$i]}{$nuc})){
					${$seq_hash{$id}}[$i]=$amb_residue;
				}
			}
		}
	}

}

#-----------------------------------------------------------------------------------

# Stretch the mask to fit a possibly gapped alignment
my $stretch_mask_ref;
my @blank_mask;
if($mask_supplied){
	if(!defined($seq_hash{$mask_seqname})){
		die "Could not find $mask_seqname in the $aln_file.\n";
	}
	my $gap_alignment_str=join "", @{$seq_hash{$mask_seqname}};
	$stretch_mask_ref=stretch_mask_to_alignment($gap_alignment_str,$mask_val_arr_ref,$gap_mask_fill_method);
	print STDERR "Mask stretched...\n";
}else{
	for(my $i=0; $i<$seq_len; $i++){
		$blank_mask[$i]=1;
	}
	$stretch_mask_ref=\@blank_mask;
}



if($dont_output_masksequenceid){
	delete($seq_hash{$mask_seqname});
	print("Removed $mask_seqname from distance matrix calculations.\n");
}

# Compute the distance matrix.  You'll want to optimize this if your matrix is really huge.
print STDERR "Computing Matrix...\n";
foreach my $a_id (sort keys %seq_hash){
	print STDERR ".";
	foreach my $b_id (sort keys %seq_hash){
		if(!defined($distance_hash{$a_id}{$b_id})){

			if(!$is_AA_alignment){
				$distance_hash{$a_id}{$b_id}=
					$distance_hash{$b_id}{$a_id}=
						nuc_score($seq_hash{$a_id}, $seq_hash{$b_id},$stretch_mask_ref);
				
			}else{
				$distance_hash{$a_id}{$b_id}=
					$distance_hash{$b_id}{$a_id}=
						aa_score($seq_hash{$a_id}, $seq_hash{$b_id},$stretch_mask_ref);

			}
		}
	}
}
print STDERR "\n";

# Dump matrix
print STDERR "Printing Matrix...\n";
printMatrix(\%distance_hash, $output_file);

print STDERR "Done.\n";

###############################################################################

sub printMatrix{
	my $matrix_hash_ref=shift;
	my $outputfilename=shift;

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

sub aa_score{
	my $a_seq_arr_ref=shift;
	my $b_seq_arr_ref=shift;
	my $stretched_mask=shift;

	my $score;

	my $seq_len_a=$#{$a_seq_arr_ref}+1;
	my $seq_len_b=$#{$b_seq_arr_ref}+1;

	if($seq_len_a!=$seq_len_b){
		die "Could not score pair, $seq_len_a != $seq_len_b\n";
	}

	#for(my $i=0; $i<$seq_len_a; $i++){
	#	print "${$a_seq_arr_ref}[$i] ${$b_seq_arr_ref}[$i] \n";
	#}

	$score=$BU->twoSequenceArrDist($a_seq_arr_ref, $b_seq_arr_ref, $stretched_mask);

	return($score);
}

#-----------------------------------------------------------------------------------

sub nuc_score{
	my $a_seq_arr_ref=shift;
	my $b_seq_arr_ref=shift;
	my $stretched_mask_ref=shift;

	my $score;

	my $seq_len_a=$#{$a_seq_arr_ref}+1;
	my $seq_len_b=$#{$b_seq_arr_ref}+1;

	if($seq_len_a!=$seq_len_b){
		die "Could not score pair, $seq_len_a != $seq_len_b\n";
	}


	my $mismatches=0;
	my $sum_mask=0;
	for(my $i=0; $i<$seq_len_a; $i++){

		my $a=${$a_seq_arr_ref}[$i];
		my $b=${$b_seq_arr_ref}[$i];
		$sum_mask+=${$stretched_mask_ref}[$i];

		if($a ne $b){
			if(($a ne "N") && ($b ne "N")){

				if(!$ignore_stutter){
					$mismatches+=${$stretched_mask_ref}[$i];	
				}else{
					# If we ignore stutter, then check both sequences 
					if(!stutter($i, $a_seq_arr_ref, $b_seq_arr_ref) && 
					   !stutter($i, $b_seq_arr_ref, $a_seq_arr_ref)){
						$mismatches+=${$stretched_mask_ref}[$i];
					}
				}
			}
		}
	}
	
	$score = $mismatches / $sum_mask;

	return($score);

}

#-----------------------------------------------------------------------------------

sub stutter{

	# Check the position on A to see if it is a gap.
	# Check the around position on B, to see if there are any nucleotides before or after the position

	my $pos=shift;
	my $a_arr_ref=shift;
	my $b_arr_ref=shift;


	# Seq A:  
	my ($nna, $na)=(${$a_arr_ref}[$pos-2], ${$a_arr_ref}[$pos-1]);
	my ($a)=${$a_arr_ref}[$pos];
	my ($an, $ann)=(${$a_arr_ref}[$pos+1], ${$a_arr_ref}[$pos+2]);

	# Seq B:
	my $b=${$b_arr_ref}[$pos];

	if($a eq "-"){
		if($b eq $nna && $b eq $na){
			return(1);
		}
		if($b eq $an && $b eq $ann){
			return(1);
		}
	}else{
		return(0);
	}

}

#-----------------------------------------------------------------------------------

sub load_mask{
	my $maskfname=shift;
	my @mask_values;

	open(FH, "<$maskfname") || die "Could not open $maskfname\n";
	while(<FH>){
		chomp;
		if(substr($_,0,1) eq "#"){
			next;
		}
		my ($res, $mask_val)=split /\t/, $_;
		push @mask_values, $mask_val;
	}
	close(FH);

	my $mask_length=$#mask_values+1;
	print STDERR "Loaded Mask Length: $mask_length\n";
	return(\@mask_values);
}

#-----------------------------------------------------------------------------------

sub print_mask_values{
	my $mask_arr_ref=shift;
	my $mask_length=$#{$mask_arr_ref}+1;
	for(my $i=0; $i<$mask_length; $i++){
		print "${$mask_arr_ref}[$i]\n";
	} 
}

#-----------------------------------------------------------------------------------

sub stretch_mask_to_alignment{
	my $alignment=shift;
	my $mask_val_arr_ref=shift;
	my $gap_fill=shift;
	my @stretched_mask;
	
	# Determine whether the gapping mask model is to fill with zeros or average the 
	# mask before and after gap
	my $ZERO=1;
	my $AVERAGE=2;
	my $gf;
	if(lc($gap_fill) eq "zero"){
		$gf=$ZERO;
	}elsif(lc($gap_fill) eq "average"){
		$gf=$AVERAGE;
	}

	my @align_arr=split //, $alignment;
	my @mask_val_arr=@{$mask_val_arr_ref};
	my $mask_len=$#mask_val_arr+1;
	my $align_len=$#align_arr+1;

	# Assert that gap free alignment length is the same as the length of the mask
	my $gap_free_alignment=$alignment;
	$gap_free_alignment=~s/-//g;
	my $gap_free_alignment_len=length($gap_free_alignment);
	if($gap_free_alignment_len!=$mask_len){
		die "Error!  Gap free alignment ($gap_free_alignment_len) doesn't match mask length ($mask_len)!\n";
	}
	
	# Go through and align mask to gapped alignment
	my $mask_pos=0;
	for(my $i=0; $i<$align_len; $i++){
		my $cur_mask_val;

		if($align_arr[$i] eq "-"){
			# If alignment has a gap
			if($gf==$ZERO){
				# Fix value to zero
				$cur_mask_val=0;

			}elsif($gf==$AVERAGE){
				my $before=$mask_pos-1;
				my $after=$mask_pos;

				# Edge cases
				if($before<0){
					$before=0;
				}
				if($after>=$mask_len){
					$after=$mask_len-1;
				}

				# Average before and after mask values
				$cur_mask_val=($mask_val_arr[$before]+$mask_val_arr[$after])/2;
			}
			push @stretched_mask, $cur_mask_val;
		}else{
			# If alignment has no gap
			push @stretched_mask, $mask_val_arr[$mask_pos];
			$mask_pos++;
		}
	}
	return(\@stretched_mask);
}

sub stretch_test{
	my @mask=(1,2,3,4,5,6);
	my $alignment="-A----BC--DEF-";
	my $stretched=stretch_mask_to_alignment($alignment, \@mask, "zero");
	print_mask_values($stretched);
	print "-------\n";
	my $stretched=stretch_mask_to_alignment($alignment, \@mask, "average");
	print_mask_values($stretched);
}
