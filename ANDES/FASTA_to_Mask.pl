#!/usr/bin/env perl

###############################################################################
#                                                                             # 
#       Copyright (c) 2011 J. Craig Venter Institute.                         #     
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
use vars qw($opt_f $opt_r $opt_o $opt_n $opt_m $opt_d);
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("f:r:o:mnd:");
my $usage = "usage: 
$0 
	-f <fasta file>
	-r <regions of interest file>
	-o <output mask name>
	[-n (normalize values all points sum to 1)]
	[-m (use region weights as max, otherwise weights are summed)]
	[-d <default value if non specified in regions file, default = 0>]

	This script will read in the fasta file and the regions of interest
	and generate a mask file.  The format of the mask file is:

		# comments
		<res>\\t<weight>\\n
		<res>\\t<weight>\\n
		...
		<res>\\t<weight>\\n

	The format of the regions of interest file should be:

		# comments
		<begin>\\t<end>\\t<weight>\\t<comments\\n
		<begin>\\t<end>\\t<weight>\\t<comments\\n
		...
		<begin>\\t<end>\\t<weight>\\t<comments\\n

	The weights can be any number.  The begin/end coordinates
	must be in 0-space based coordinates.
	
	Using the -n flag will force all the values to equal 1 when summed.

	Using the -m flag will take the maximum weight if there are overlapping
	regions.  By default, the sum of all the weights is used.

	-n and -m are not mutually exclusive.

	Note if you set a default value with the -d, the normalization will
	change the default value to be relative to your other values.

";

if(!defined($opt_f) || !defined($opt_r) || !defined($opt_o)){
	die $usage;
}

my $FastaFile=$opt_f;
my $RegionsFile=$opt_r;
my $OutputMaskFile=$opt_o;

my $DoNormalize=0;
if(defined($opt_n)){
	$DoNormalize=1;
}

my $DoMax=0;
if(defined($opt_m)){
	$DoMax=1;
}

my $DefaultVal=0;
if(defined($opt_d)){
	$DefaultVal=$opt_d;
}


###############################################################################

print STDERR "Input FASTA: $FastaFile\n";
print STDERR "Regions File: $RegionsFile\n";
print STDERR "Output Mask: $OutputMaskFile\n";
print STDERR "Normalize?: $DoNormalize\n";
print STDERR "Take Max?: $DoMax\n";
print STDERR "Default Value: $DefaultVal\n";
print STDERR "\n";

###############################################################################
# Load regions function

sub load_regions{
	my $fname=shift;
	my @begins;
	my @ends;
	my @weights;

	open(FH, "<$fname") || die "Could not open $fname\n";
	while(<FH>){
		chomp;
		if(substr($_,0,1) eq "#"){
			next;
		}
		my ($begin, $end, $weight, $other)=split /\t/, $_;
		push @begins, $begin;
		push @ends, $end;
		push @weights, $weight;
	}
	close(FH);
	return(($#begins+1), \@begins, \@ends, \@weights);
}

###############################################################################
# Load sequence from fasta file

my @ref_seq_arr;

open(FH, "<$FastaFile") || die "Could not open $FastaFile\n";

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
		$_=~s/\r//g;
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

close(FH);

#------------------------------------------------------------------------------

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	@ref_seq_arr=split //, $sequence;
}

###############################################################################

my ($num_regions, $begins_ref, $ends_ref, $weights_ref)=load_regions($RegionsFile);
print STDERR "Num Regions loaded: $num_regions\n";

# Set up empty mask value array
my @mask_values;
my @touched;
my $seq_len=$#ref_seq_arr+1;
print STDERR "Sequence length: $seq_len\n";
# Initialize mask values
for(my $i=0; $i<$seq_len; $i++){
	$mask_values[$i]=0;
	$touched[$i]=0;
}


# Fill in mask with region values
for(my $regix=0; $regix<$num_regions; $regix++){
	my $reg_begin= @{$begins_ref}[$regix];		
	my $reg_end=   @{$ends_ref}[$regix];		
	my $reg_weight=@{$weights_ref}[$regix];		

	for(my $i=$reg_begin; $i<$reg_end; $i++){
		if($DoMax){
			# Take the mask value if regions overlap
			if($mask_values[$i]<$reg_weight){
				$mask_values[$i]=$reg_weight;
			}
		}else{
			# Sum up mask values if regions overlap
			$mask_values[$i]+=$reg_weight;
		}
		$touched[$i]=1;
	}
}

# Set values for regions that have not be specified
for(my $i=0; $i<$seq_len; $i++){
	if(!$touched[$i]){
		$mask_values[$i]=$DefaultVal;
	}
}

# Perform normalization if requested
if($DoNormalize){
	print STDERR "Performing normalization.\n";

	# Sum to compute normalization factor
	my $sum=0;
	for(my $i=0; $i<$seq_len; $i++){
		$sum+=$mask_values[$i];
	}

	# Perform normalization
	for(my $i=0; $i<$seq_len; $i++){
		$mask_values[$i]/=$sum;
	}
}

# Output
print STDERR "Writing mask file out.\n";
open(FH, ">$OutputMaskFile") || die "Could not open $OutputMaskFile\n";
for(my $i=0; $i<$seq_len; $i++){
	print FH "$ref_seq_arr[$i]\t$mask_values[$i]\n";	
}
close(FH);

#------------------------------------------------------------------------------

print STDERR "Completed.\n";
