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
use vars qw($opt_a $opt_p);
use File::Basename;
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;
use BLOSUM_Utilities;

getopts("a:p:");
my $usage = "usage: 
$0 
	-a <.ALN file>
	-p <positions file>

	Takes aln file from clustalw and a list of positions and extracts the column/alleles
	and places them into a table.  

	The positions file positions should start from 1 and go to N.
	
";

if(!defined($opt_a) || !defined($opt_p)){
	die $usage;
}

my $aln_file=$opt_a;
my $pos_file=$opt_p;

my $position_output_file=$pos_file;
$position_output_file=~s/txt$/alleles/;

###############################################################################

my $seq_hash_ref=PositionProfiles::readClustalAlnFile($aln_file);

###############################################################################

my @positions;
open(FH, "<$opt_p") || die "Could not open $opt_p\n";
while(<FH>){
	chomp;
	push @positions, $_;
}
close(FH);
my $num_positions=$#positions+1;

print STDERR "Num Positions: $num_positions\n";

###############################################################################
# Output alleles

my @samples=keys %{$seq_hash_ref};
my $num_samples=$#samples+1;

print STDERR "Num Samples: $num_samples\n";

open(ALLELE_FH, ">$pos_file\.alleles") || die "Could not open $pos_file\.alleles\n";
print ALLELE_FH "sequence_id\t" . (join "\t", @positions) . "\n";
for(my $i=0; $i<$num_samples; $i++){
	my @alignment=@{${$seq_hash_ref}{$samples[$i]}};
	
	my @alelles;
	for(my $j=0; $j<$num_positions; $j++){
		push @alelles, $alignment[$positions[$j]-1];
	}

	print ALLELE_FH "$samples[$i]\t";
	print ALLELE_FH join "\t", @alelles;
	print ALLELE_FH "\n";
}
close(ALLELE_FH);

###############################################################################
# Output distance from median
my $peptide_similarity_matrix="$FindBin::Bin/BLOSUM/BLOSUM62";
my $BU=new BLOSUM_Utilities($peptide_similarity_matrix);

my @positional_distances;

# Compute distance between most common residue versus other residues in position
for(my $i=0; $i<$num_positions; $i++){

	my $pos=$positions[$i];
	print STDERR "Position: $pos\n";

	# Compute residue frequencies
	my %pos_hash;
	for(my $j=0; $j<$num_samples; $j++){
		my $sample_id=$samples[$j];
		my @alignment_arr=@{${$seq_hash_ref}{$sample_id}};
		my $residue=$alignment_arr[$pos-1];
		#print STDERR "Residue: $residue\n";
		$pos_hash{$residue}++;
	}	

	# Find mode/most common residue
	my $most_common_res="";
	my $most_common_freq=0;
	foreach my $res(keys %pos_hash){
		if($pos_hash{$res}>$most_common_freq){
			$most_common_freq=$pos_hash{$res};
			$most_common_res=$res;
		}
		print "  $res: $pos_hash{$res}\n";
	}
	print STDERR "  Mode Residue: $most_common_res\n";

	# Compute distance to mode residue
        for(my $samp_idx=0; $samp_idx<$num_samples; $samp_idx++){
		my $sample_id=$samples[$samp_idx];
		my @alignment_arr=@{${$seq_hash_ref}{$sample_id}};
		my $residue=$alignment_arr[$pos-1];
		my $score=1-$BU->similarity_score($most_common_res, $residue);
		#print STDERR "$residue: $score\n";
		$positional_distances[$samp_idx][$i]=$score;
	}

	print STDERR "\n";
}

# Output distances
open(DIST_FH, ">$pos_file\.dist") || die "Could not open $pos_file\.dist\n";
print DIST_FH "sequence_id\t" . (join "\t", @positions) . "\n";
for(my $samp_idx=0; $samp_idx<$num_samples; $samp_idx++){
	print DIST_FH "$samples[$samp_idx]";
	for(my $pos_idx=0; $pos_idx<$num_positions; $pos_idx++){
		printf DIST_FH "\t%2.4f", $positional_distances[$samp_idx][$pos_idx];
	}
	print DIST_FH "\n";
}
close(DIST_FH);

###############################################################################

print STDERR "Done!\n";

###############################################################################
