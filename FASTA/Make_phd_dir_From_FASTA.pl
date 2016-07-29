#!/usr/bin/env perl

###############################################################################

use FindBin ();
use lib "$FindBin::Bin";

use strict;
use Getopt::Std;
use vars qw($opt_s $opt_q $opt_o);

getopts("s:q:o:");
my $usage = "usage:
$0
        -s <sequence fasta file>
        -q <quality fasta file>
	-o <output dir>

	This program will take in a sequence and quality multi fasta file
	and generate the files necessary for consed to work properly.  Consed
	is very picky about how files are generated, so there are a lot of
	ridiculous steps that seems to be inefficient in this script.

	Once this script has completed, you should be able to do something like:

	cap3 <output dir>/asm_dir/phd_based.fasta

	and then the output from cap3: <output dir>/asm_dir/phd_based.fasta.cap.ace

	should be viewable with consed with the quality values properly displayed.
	
";

if(!defined($opt_s) || !defined($opt_q) || !defined($opt_o)){
        die $usage;
}

my $seq_fn=$opt_s;
my $qual_fn=$opt_q;
my $output_dir=$opt_o;

print STDERR "Making necessary directories in: $output_dir\n";
if(!(-e $output_dir)){
	mkdir $output_dir;
	mkdir  "$output_dir/phd_dir";
	mkdir  "$output_dir/chromat_dir";
	mkdir "$output_dir/asm_dir";
}else{
	die "This program shouldn't be run in an existing directory.  It could be sensitive and destructive" .
		" to what's currently in there.\n";
}

print STDERR "Splitting input fasta files into individual single record fastas.\n";
`$FindBin::Bin/Split_MultiFASTA_into_FASTAs.pl -f $seq_fn -d $output_dir/seq`;
`$FindBin::Bin/Split_MultiFASTA_into_FASTAs.pl -f $qual_fn -d $output_dir/qual`;

print STDERR "Running 'make trace':\n";
foreach my $seq(split /\n/, `ls -1 $output_dir/seq`){
	print STDERR "\t$seq\n";
	my $changed_name=$seq;
	$changed_name=~s/\.fasta$/.seq/;
	`mv $output_dir/qual/$seq $output_dir/$changed_name.qual`;
	`mv $output_dir/seq/$seq $output_dir/$changed_name`;
	`cd $output_dir; mktrace $changed_name $changed_name\.scf`;
}

print STDERR "Moving scf and phd files to appropriate chromat and phd directories.\n";
`mv $output_dir/*.scf $output_dir/chromat_dir`;
`mv $output_dir/*.phd.1 $output_dir/phd_dir`;

print STDERR "Running phd2fasta to generate assembly input files with proper deflines.\n";
`cd $output_dir; phd2fasta -id phd_dir -os phd_based.fasta -oq phd_based.fasta.qual`;

print STDERR "Moving assembly input files to asm_dir.\n";
`mv $output_dir/phd_based.fasta $output_dir/phd_based.fasta.qual $output_dir/asm_dir`;

print STDERR "Removing temporary seq and qual directories and intermediate single record fastas.\n";
rmdir "$output_dir/qual/";
rmdir "$output_dir/seq/";
`rm $output_dir/*.seq $output_dir/*.qual`;

print STDERR "Done.\n";



