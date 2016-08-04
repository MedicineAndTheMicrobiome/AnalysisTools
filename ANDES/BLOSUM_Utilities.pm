#!/usr/bin/perl

###############################################################################
#                                                                             #
#       Copyright (c) 2005 J. Craig Venter Institute.                         #
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

package BLOSUM_Utilities;

use Carp;
use strict;

my %blosum_hash;
my %similarity_hash;

###############################################################################

sub new {
	my $this = shift;
	my $matrix_name = shift;
	my $gap_penalty = shift;

	my $class = ref($this) || $this;
	my $self = {};
	bless $self, $class;

	if(!defined($gap_penalty)){
		$self->{gap_penalty}=-10;		
	}

	$self->loadBLOSUM($matrix_name);
	$self->compute_similarity_hash();

	return $self;
}

###############################################################################

sub loadBLOSUM{
	my $self = shift;
	my $filename=shift;

	print STDERR "Loading BLOSUM Matrix: $filename\n";
	open(FH, "<$filename") || die "Could not open $filename.\n";

	my @header;
	while(<FH>){
		chomp;
		if($_=~/^#/){
			next;
		}else{
			
			if($_=~/^ /){
				$_=~s/^\s+//;
				@header=split /\s+/, $_;		
			}else{
				my @columns=split /\s+/, $_;
				for(my $i=1; $i<=$#columns; $i++){
					$blosum_hash{$columns[0]}{$header[$i-1]}=@columns[$i];
				}
			}		
		}
	}

	close(FH);

	# Insert gap penalty
	foreach my $res(keys %blosum_hash){
		$blosum_hash{"-"}{$res}=$self->{gap_penalty};
		$blosum_hash{$res}{"-"}=$self->{gap_penalty};
	}
	$blosum_hash{"-"}{"-"}=1;
	
	print STDERR "done.\n";
}

###############################################################################

sub compute_similarity_hash{
	my $self = shift;

	my @residues=keys %blosum_hash;

	foreach my $res1(@residues){
		my $min=0;
		my $max=0;
		my $range;

		# Find min/max/range
		foreach my $res2(@residues){
			my $score=$blosum_hash{$res1}{$res2};

			if($min>$score){
				$min=$score;
			}elsif($max<$score){
				$max=$score;
			}

			$range=$max-$min;			
		}

		# Normalize BLOSUM matrix into similarity matrix
		foreach my $res2(@residues){
			my $score=$blosum_hash{$res1}{$res2};
			my $adj=($score-$min)/$range;
			$similarity_hash{$res1}{$res2}=$adj;
		}
		
	}
	
}

###############################################################################

sub print_matrix{
	my $self=shift;

	my @residues=sort(keys %blosum_hash);

	print "   " . (join "   ", ("",@residues));
	print "\n";
	for(my $i=0; $i<=$#residues; $i++){
		print("$residues[$i]  ");
		for(my $j=0; $j<=$#residues; $j++){
			printf("%4i", $blosum_hash{$residues[$i]}{$residues[$j]});
		}
		print "\n";
	}

	print "\n\n";

	print "  " . (join "    ", ("",@residues));
	print "\n";
	for(my $i=0; $i<=$#residues; $i++){
		print("$residues[$i]  ");
		for(my $j=0; $j<=$#residues; $j++){
			printf("%3.2f ", $similarity_hash{$residues[$i]}{$residues[$j]});
		}
		print "\n";
	}

}

###############################################################################

sub blosum_score{
	my $self=shift;
	my $a=shift;
	my $b=shift;

	if(defined($blosum_hash{$a}{$b})){
		return($blosum_hash{$a}{$b});
	}else{
		die "Undefined AA comparison: '$a' <=> '$b'\n";
	}
}

sub similarity_score{
	my $self=shift;
	my $a=shift;
	my $b=shift;

	
	if(defined($similarity_hash{$a}{$b})){
		return($similarity_hash{$a}{$b});
	}else{
		die "Undefined AA comparison: '$a' <=> '$b'\n";
	}
}
	
###############################################################################

sub twoSequenceArrScore{
	my $self = shift;
	my $seq_arr1_ref = shift;
	my $seq_arr2_ref = shift;
	my $mask_ref = shift;

	my $mask_specified=0;
	my @mask_values;
	if(defined($mask_ref)){
		print STDERR "Mask specified.\n";
		$mask_specified=1;
		@mask_values=@{$mask_ref};
	}

	my $seq1len=$#{$seq_arr1_ref}+1;
	my $seq2len=$#{$seq_arr2_ref}+1;
	if($seq1len != $seq2len){
		die "Cannot compute two sequence distance.  Lengths do not match: $seq1len != $seq2len . \n";
	}

	if($mask_specified){
		my $masklen=$#mask_values+1;
		if($masklen!=$seq1len){
			die "Error:  Mask length does not match (aligned) sequence lengths. Mask: $masklen  Seq: $seq1len\n";
		}
	}else{
		@mask_values=split //, ("1" x $seq1len);
	}

	my $total_blosum_score=0;
	my $compare_length=0;
	for(my $i=0; $i<$seq2len; $i++){
		if((${$seq_arr1_ref}[$i] ne " ") && (${$seq_arr2_ref}[$i] ne " ")){
			$total_blosum_score+=$self->blosum_score(${$seq_arr1_ref}[$i], ${$seq_arr2_ref}[$i])*$mask_values[$i];

			if($mask_values[$i]>0){
				$compare_length++;
			}
		}
	}

	return($total_blosum_score);
	
}

sub twoSequenceScore{
	my $self = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $mask_ref = shift;

	my @seq1arr=split //, $seq1;
	my @seq2arr=split //, $seq2;

	return($self->twoSequenceArrScore(\@seq1arr, \@seq2arr, $mask_ref));

}

###############################################################################

sub twoSequenceArrDist{
	my $self = shift;
	my $seq_arr1_ref = shift;
	my $seq_arr2_ref = shift;
	my $mask_ref = shift;

	my $mask_specified=0;
	my @mask_values;
	if(defined($mask_ref)){
		$mask_specified=1;
		@mask_values=@{$mask_ref};
	}

	my $seq1len=$#{$seq_arr1_ref}+1;
	my $seq2len=$#{$seq_arr2_ref}+1;
	if($seq1len != $seq2len){
		die "Cannot compute two sequence distance.  Lengths do not match: $seq1len != $seq2len . \n";
	}

	if($mask_specified){
		my $masklen=$#mask_values+1;
		if($masklen!=$seq1len){
			die "Error:  Mask length does not match (aligned) sequence lengths. Mask: $masklen  Seq: $seq1len\n";
		}
	}else{
		@mask_values=split //, ("1" x $seq1len);
	}

	#my $str1=join "", @{$seq_arr1_ref};
	#my $str2=join "", @{$seq_arr2_ref};
	#print STDERR "1: '$str1'\n";
	#print STDERR "2: '$str2'\n";
	#print STDERR "\n\n";

	my $total_similarity=0;
	my $compare_length=0;
	for(my $i=0; $i<$seq2len; $i++){
		#print STDERR "${$seq_arr1_ref}[$i] / ${$seq_arr2_ref}[$i]\n";
		if((${$seq_arr1_ref}[$i] ne " ") && (${$seq_arr2_ref}[$i] ne " ")){
			$total_similarity+=$self->similarity_score(${$seq_arr1_ref}[$i], ${$seq_arr2_ref}[$i])*$mask_values[$i];
			$compare_length+=$mask_values[$i];
		}
	}

	#print STDERR "$compare_length\n";
	if($compare_length > 0){
		return(1-($total_similarity/$compare_length))
	}else{
		return(1);
	}
	
}

sub twoSequenceDist{
	my $self = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $mask_ref = shift;

	my @seq1arr=split //, $seq1;
	my @seq2arr=split //, $seq2;

	return($self->twoSequenceArrDist(\@seq1arr, \@seq2arr, $mask_ref));

}


###############################################################################

1;
