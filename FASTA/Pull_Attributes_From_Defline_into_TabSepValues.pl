#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;

###############################################################################

use vars qw ($opt_f $opt_t);
getopts("t:f:p:");

my $usage = "
	$0
	    -f <FASTA file>
	    -t <comma separated list of tags in fasta file>

	for example:

	$0 -f myfasta -t clear_start,clear_end

	would produce:

	1117995667794\\t10\\t20\\n

	from a fasta with a defline that looks like:

	>1117995667794 /clear_start=10 /clear_end=20
	
";

if(!defined($opt_f) || !defined($opt_t)){
	die $usage;
}

my $fasta_fn=$opt_f;
my @taglist=split /,/, $opt_t;

###############################################################################

my %attributes_hash;

###############################################################################

my ($defline, $prev_defline, $sequence);
open(FASTA_FH, "<$fasta_fn") || die "Could not open $fasta_fn\n";
print STDERR "Processing FASTA file...\n";
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
                $sequence.=$_;
        }
}
process_record($prev_defline, $sequence);

print STDERR "\n";
print STDERR "done.\n";

close(FASTA_FH);

sub process_record{
	my $defline=shift;
	my $sequence=shift;

	my $id;
	my $len=length($sequence);

	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse out ID, $_\n"
	}

	my @attributes=($id);
	foreach my $tag(@taglist){
		if($defline=~/$tag=(\S+)/){
			push @attributes,$1;
		}else{
			push @attributes,"";
		}
	}

	print STDOUT (join "\t", @attributes) . "\n";
	
}

