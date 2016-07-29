#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	-f <fasta file>

	This program will read in a multi-fasta file and sum up the total
	number of residues in the file.
";

if(!(
	defined($opt_f))){
	die $usage;
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

###############################################################################

print STDERR "Reading in FASTA file...\n";
my ($defline, $prev_defline, $sequence);
my $i;

my $sum_seq_lengths=0;

$sequence="";
while(<FASTA_FH>){
        chomp;
        if(/^>/){
                $defline=$_;
                if($sequence ne ""){
			$sum_seq_lengths+=process_record($prev_defline, $sequence);
                        $sequence="";
                }
                $prev_defline=$defline;
        }else{
                $sequence.="$_";
        }
}
$sum_seq_lengths+=process_record($prev_defline, $sequence);

close(FASTA_FH);

print "$sum_seq_lengths\n";
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline=shift;
	my $sequence=shift;
	
	return(length($sequence));
	
}


