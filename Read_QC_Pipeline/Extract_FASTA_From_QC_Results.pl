#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_f $opt_r $opt_p $opt_c);
use File::Basename;

my $fastq_splitter_bin="$FindBin::Bin/../FASTQ/Split_Fastq/Split_Fastq.pl";

getopts("i:o:frpc");
my $usage = "usage: 

$0 

	-i <Input Directory, i.e. Output directory of QC Pipeline>
	-o <Output directory>

	-f (Extract forward)
	-r (Extract reverse)
	-p (Extract paired)
	-c (Combine Reads into Single File)

	This script will take the output processed FASTQ files
	and extract them into FASTA files, saving them into the
	output directory.

	The -f, -r, and -p option specify which reads to extract.
	
	If the -c option is specified, all the reads will be combined
	into a single output FASTA file.

	This script will use the FASTQ utility:
	$fastq_splitter_bin
	
";

if(!defined($opt_i)||!defined($opt_o)){
	die $usage;
}

my $input_dir=$opt_i;
my $output_dir=$opt_o;

my $extr_for;
my $extr_rev;
my $extr_pair;
my $combine;

$extr_for=defined($opt_f)?1:0;	
$extr_rev=defined($opt_r)?1:0;	
$extr_pair=defined($opt_p)?1:0;	
$combine=defined($opt_c)?1:0;	

###############################################################################

print STDERR "\n";
print STDERR "Input Dir: $input_dir\n";
print STDERR "Output Dir: $output_dir\n";
print STDERR "\n";
print STDERR "Extract:\n";
print STDERR "  Forward: $extr_for\n";
print STDERR "  Reverse: $extr_rev\n";
print STDERR "  Paired: $extr_pair\n";
print STDERR "\n";
print STDERR "Combine?: $combine\n";
print STDERR "\n";

###############################################################################

if(!(-e $output_dir)){
	print STDERR "Making Output Directory: $output_dir\n";
	mkdir $output_dir;
}
if(!(-e $output_dir)){
	die "Error Making $output_dir\n";
}

###############################################################################

print STDERR "Getting Sample Directory List...\n";

my $cwd=`pwd`;

my $list=`cd $input_dir; ls -d */; cd $cwd`;

$list=~s/\/\n/\n/g;

my @samp_dir_list=split /\n/, $list;

foreach my $sample (@samp_dir_list){
	print STDERR "$sample\n";

	my $for_fn="$input_dir/$sample/$sample.for.frag.fastq";
	my $rev_fn="$input_dir/$sample/$sample.rev.frag.fastq";
	my $pair_fn="$input_dir/$sample/$sample.paired.fastq";

	if(!$combine){
		if($extr_for){
			`$fastq_splitter_bin -f $for_fn -a $output_dir/$sample.for.frag.fasta`;
		}
		if($extr_rev){
			`$fastq_splitter_bin -f $rev_fn -a $output_dir/$sample.rev.frag.fasta`;
		}
		if($extr_pair){
			`$fastq_splitter_bin -f $pair_fn -a $output_dir/$sample.paired.fasta`;
		}
	}else{

		`touch $output_dir/$sample.fasta`;

		if($extr_for){
			`$fastq_splitter_bin -f $for_fn -a $output_dir/tmp`;
			`cat $output_dir/tmp >> $output_dir/$sample.fasta`;
		}
		if($extr_rev){
			`$fastq_splitter_bin -f $rev_fn -a $output_dir/tmp`;
			`cat $output_dir/tmp >> $output_dir/$sample.fasta`;
		}
		if($extr_pair){
			`$fastq_splitter_bin -f $pair_fn -a $output_dir/tmp`;
			`cat $output_dir/tmp >> $output_dir/$sample.fasta`;
		}

		`rm $output_dir/tmp`;
	}
}


###############################################################################


sub convert_fastq_to_fasta{
	my $fastq_path=shift;
	my $config=shift;

	# Generate output fasta file name
	my $fasta_path=$fastq_path;
	$fasta_path=~s/\.fastq.gz$/.fasta/;
	$fasta_path=~s/\.fastq$/.fasta/;
	$fasta_path=~s/\.fq.gz$/.fasta/;
	$fasta_path=~s/\.fq$/.fasta/;

	my $path=$config->val("FASTQ_Tools", "Path");
	my $bin=$config->val("FASTQ_Tools", "SplitFASTQ_bin");

	print STDERR "Running Convert FASTQ to FASTA: $bin\n";
	my $res=`$path/$bin -f $fastq_path -a $fasta_path`;

	return($fasta_path);
}


###############################################################################

print STDERR "Done.\n";

