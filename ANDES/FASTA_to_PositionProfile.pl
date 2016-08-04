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
use vars qw($opt_f $opt_o $opt_m);
use FindBin ();
use lib $FindBin::Bin;
use PositionProfiles;

getopts("f:o:m");
my $usage = "usage: 
$0 
	-f <fasta file>
	[-o <output profile name>]
	[-m (aMino acid flag)]

	Writes out a profile for a single fasta.

	If the input fasta file has ambiguity codes it in, the underlying
	bases will be fairly represented.

	This program is NOT for generating proper profiles based on a multifasta
	of similar sequences.

	This program may be used for creating a reference profile based on a 
	reference fasta.  Or for generating a profile for a single fasta (that
	could not be clustalw'd).

";

if(!defined($opt_f)){
	die $usage;
}

my $outfile;
if(defined($opt_o)){
	$outfile=$opt_o;
}

my $residue_type=PositionProfiles::NUC_TYPE;
if(defined($opt_m)){
	$residue_type=PositionProfiles::AA_TYPE;
}

print STDERR "Residue Type: $residue_type\n";

my $fastafn=$opt_f;

###############################################################################

open(FH, "<$fastafn") || die "Could not open $fastafn\n";

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
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my @pos_arr;
	
	my @nucs=split //, uc($sequence);
	for(my $i=0; $i<=$#nucs; $i++){
		PositionProfiles::ambig_to_prof($nucs[$i], \%{$pos_arr[$i]}, $residue_type);	
	}

	my $filename;
	if(defined($outfile)){
		$filename=$outfile;
	}else{
		if($defline=~/^>(\S+)/){
			$filename="$1\.$residue_type\.prof";
		}else{
			die "Error parsing defline for output profile name.\n";
		}
	}

	PositionProfiles::output_prof($filename, \@pos_arr, $sequence, $residue_type);
}

