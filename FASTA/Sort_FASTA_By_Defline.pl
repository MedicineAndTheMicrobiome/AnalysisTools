#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	-f <fasta file>

	This program will read in the specified FASTA file and then sort
	the records by their deflines and then send a sorted FASTA file
	to STDOUT.  The width of each output record sequence will stay the
	same.

	This is a memory hog because the entire FASTA file is stored in a 
	hash.
";

if(!(
	defined($opt_f))){
	die $usage;
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

###############################################################################

my @defline_array;
my %sequence_hash;
my $i=0;

print STDERR "Reading in FASTA file...\n";
my ($defline, $prev_defline, $sequence);
$sequence="";
while(<FASTA_FH>){
        chomp;
        if(/^>/){
                $defline=$_;
                if($sequence ne ""){
                        process_record($prev_defline, $sequence, $i);
                        $sequence="";
			$i++;
                }
                $prev_defline=$defline;
        }else{
                $sequence.="$_\n";
        }
}
process_record($prev_defline, $sequence, $i);

close(FASTA_FH);

# Sort defline array
print STDERR "Sorting by defline...\n";
@defline_array=sort sort_function @defline_array;

print STDERR "Outputing results...\n";
for($i=0;$i<=$#defline_array;$i++){
	print "$defline_array[$i]\n";
	$sequence=$sequence_hash{$defline_array[$i]};
	print $sequence;
}

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline=shift;
	my $sequence=shift;
	my $i=shift;
	
	$defline_array[$i]=$defline;
	$sequence_hash{$defline}=$sequence;

}

###############################################################################

sub sort_function{
	$a cmp $b;
}

###############################################################################

sub grab_from_defline{
        my $defline=shift;
        my $key=shift;
		my $i;

		my @tags=split /\/(\S+)=/, $defline;
		for($i=1; $i<$#tags; $i+=2){
			if($key eq $tags[$i]){
				return($tags[$i+1]);
			}
		}
        return("");
}


