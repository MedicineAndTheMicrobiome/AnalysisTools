#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_d);

getopts("f:d:");
my $usage = "usage: 
$0 
	-f <fasta file>
	[-d <output dir>]

	This program will read in a multi-FASTA file and generate multiple FASTA files.
	Output name will be based on the defline.

	This program is safe for quality files.
";

if(!(
	defined($opt_f))){
	die $usage;
}

my $out_dir=".";
if(defined($opt_d)){
	$out_dir=$opt_d;
	if(!(-e $out_dir)){
		mkdir $out_dir;
	}
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

###############################################################################

print STDERR "Reading in FASTA file...\n";
my ($defline, $prev_defline, $sequence);
my $i;
$sequence="";
while(<FASTA_FH>){
        chomp;
        if(/^>/){
                $defline=$_;
                if($sequence ne ""){
                        process_record($prev_defline, $sequence);
                        $sequence="";
                }
                $prev_defline=$defline;
        }else{
                $sequence.="$_\n";
        }
}
process_record($prev_defline, $sequence);

close(FASTA_FH);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline=shift;
	my $sequence=shift;
	
	if($defline=~/^>(\S+)/){
		my $id = $1 . ".fasta";
		open(OUT_FH, ">$out_dir/$id") || die "Could not open $out_dir/$id for writing...\n";
		print OUT_FH "$defline\n";
		print OUT_FH "$sequence";
		close(OUT_FH);
	}else{
		die "Could not parse defline for filename: $defline\n";
	}
}


