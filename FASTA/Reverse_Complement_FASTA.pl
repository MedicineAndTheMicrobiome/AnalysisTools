#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_l);

getopts("l:");

my $usage = "usage: 
$0 
	[-l <ids in fasta stream to reverse complement, optional>]

	This program will read in a FASTA file through STDIN and then
	output the same FASTA file through STDOUT, except with the 
	sequences reverse complemented.  If the -l option is specified
	then only the ids in the will will be reverse complemented.
";

print STDERR "$usage\n";

my $rc_list=$opt_l;
my $use_rc_list=defined($rc_list);

###############################################################################

my %rc_id_hash;

if($use_rc_list){

	open(LIST_FH, "<$rc_list") || die "Could not open $rc_list for reading.\n";

	while(<LIST_FH>){
		chomp;
		my ($id, $otherstuff)=split /\t/, $_;
		$rc_id_hash{$id}=1;
	}

	close(LIST_FH);
}

###############################################################################

print STDERR "Processing FASTA file...\n";

my $rc_count=0;
my $total_count=0;
my ($defline, $prev_defline, $sequence);
while(<STDIN>){
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

print STDERR "$rc_count of $total_count sequences have been reverse complemented.\n";
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse defline for sequence id: $defline\n";
	}

	$total_count++;
	if(!$use_rc_list || ($use_rc_list && defined($rc_id_hash{$id}))){
		$sequence=revcomp($sequence);
		$rc_count++;
	}
	
	#print STDERR "Working on $defline\n";
	print STDOUT "$defline\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------

sub revcomp{
    my $sequence = shift;

    $sequence =~ tr/atugcyrswkmbdhvnxATUGCYRSWKMBDHVNX/taacgryswmkvhdbnxTAACGRYSWMKVHDBNX/;
    my $revcomp = "";
    for (my $j=length($sequence)-1;$j>-1;$j--) {
        $revcomp .= substr($sequence,$j,1);
    }
    return $revcomp;

}

