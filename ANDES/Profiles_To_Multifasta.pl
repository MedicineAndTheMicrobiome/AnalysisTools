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
use vars qw($opt_o);
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use PositionProfiles;

getopts("o:");
my $usage = "usage: 
$0 
	-o <output multifasta file> <profile1> ... <profilen>

	Reads in n profile files and generates a multifasta file with n
	records.  Each sequence record in the multifasta file will be based on
	the major allele that was specified in the position profile.

	The sequence id used in the fasta defline will be based on the
	name of the profile.
	
";

if(!defined($opt_o)){
	die $usage;
}

my $output_file=$opt_o;

open(OUT_FH, ">$output_file") || die "Could not open output fasta: $output_file\n";

my $prof_name;
while($prof_name=shift){

	my ($name)=fileparse($prof_name);

	my ($major_res_arr_ref, $prof_arr_ref)=PositionProfiles::load_prof($prof_name);
	
	my $major_seq_str=join "", @{$major_res_arr_ref};
	$name=~s/\.prof$//;
	output_fasta($name, $major_seq_str);

}

close(FH);

print STDERR "Done.\n";


###############################################################################

sub output_fasta{
	my $seqname=shift;
	my $sequence=shift;

	print OUT_FH ">$seqname\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print OUT_FH substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

}
