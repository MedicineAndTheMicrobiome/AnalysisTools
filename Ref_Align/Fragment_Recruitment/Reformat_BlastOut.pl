#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_b $opt_o);
use FileHandle;

getopts("b:o:");

my $usage = "
	usage:
	$0
		-b <input blast output>
		[-o <output root>]

	
	This script will take the output from a \"-m 8\" tab delimited blast output
	file and extract the columns necessary for the Plot_FR.r script.

	If there was more than one subject in the blast file, then more
	than one output file will be generated.

";

if(!defined($opt_b)){
	die $usage;
}

my $blastfile=$opt_b;
my $output_root=$blastfile;

$output_root=~s/\.blast$//;
$output_root=~s/\.out$//;

if(defined($opt_o)){
	$output_root=$opt_o;
}

print "Blast output in file: $blastfile\n";
print "Output Root: $output_root\n";

###############################################################################

open(IN_FH, "<$blastfile") || die "Could not open map file $blastfile\n";

my %fh_hash;

my $BEGIN=9-1;
my $END=10-1;
my $PERC_ID=3-1;
my $SUBJ_ID=2-1;

while(<IN_FH>){
	chomp;
	my @in=split /\t/, $_;
	
	my ($subj_id, $begin, $end, $perc_id)=($in[$SUBJ_ID], $in[$BEGIN], $in[$END], $in[$PERC_ID]);
	
	if(!defined($fh_hash{$subj_id})){
		my $new_fname="$output_root\.$subj_id\.pos";
		print STDERR "Opening file for writing: $new_fname\n";
		$fh_hash{$subj_id}=FileHandle->new($new_fname, "w");
	}

	print {$fh_hash{$subj_id}} "$begin\t$end\t$perc_id\n";

}

close(IN_FH);

foreach my $fh(keys %fh_hash){
	print STDERR "Closing file: $fh\n";
	close($fh);
}

