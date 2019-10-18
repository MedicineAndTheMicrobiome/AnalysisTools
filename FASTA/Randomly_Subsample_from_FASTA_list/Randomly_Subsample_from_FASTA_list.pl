#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_l $opt_n $opt_r);

getopts("l:n:r:");
my $usage = "usage: 
$0 
	-l <fasta list file name>
	-n <sample size, n>
	-r <output directory>

	This program will randomly sample n sequences from each
	fasta file specified in the list.  If the number of 
	records is less than n, then just return the original
	fasta file.  Replacement will not be performed to get
	the requested sample size.
	
	The format of the list should be two columns:

	Two columns:
	<sample name>\\t<fasta filename>\\n

";

if(!(
	defined($opt_l) && 
	defined($opt_n) && 
	defined($opt_r))){
	die $usage;
}

###############################################################################

my $fasta_list=$opt_l;
my $sample_size=$opt_n;
my $output_root=$opt_r;

print STDERR "\n";
print STDERR "Input fasta list: $fasta_list\n";
print STDERR "Sample size: $sample_size\n";
print STDERR "Output root: $output_root\n";
print STDERR "\n";

if(!(-e $output_root)){
	print STDERR "Making $output_root\n";
	mkdir $output_root;
}

open(LST_FH, "<$fasta_list") || die "Could not open $fasta_list\n";

print STDERR "Loading FASTA list...\n";

my %sample_to_fasta_hash;
while(<LST_FH>){
	chomp;
	my @fields=split /\t/, $_;
	$sample_to_fasta_hash{$fields[0]}=$fields[1];
}	
		
close(LST_FH);


foreach my $sampid (keys %sample_to_fasta_hash){

	my $fasta_fn=$sample_to_fasta_hash{$sampid};	
	print STDERR "Sample ID: $sampid, File Name: $fasta_fn\n";

	# Determine num reads in fasta file
	print STDERR "Counting number of records in FASTA file\n";
	my $num_records=count_records($fasta_fn);
	print STDERR "Number of records: $num_records\n";

	# Determine if we need to subsample or not
	my $select_arr=undef;
	if($num_records<=$sample_size){
		print STDERR "Requested sample size less than or equal to number of available sequences\n";
	}else{
		# Uniformly sample from sample_size without replacement
		$select_arr=select_samples($sample_size, $num_records);
	}

	# Open input fasta file
	my $infh;
	print STDERR "Reading: $fasta_fn\n";
	open($infh, "<$fasta_fn") || die "Could not open $fasta_fn\n";

	# Open output fasta file
	my $outfn="$output_root/$sampid\.fasta";
	my $outfh;
	print STDERR "Writing: $outfn\n";
	open($outfh, ">$outfn") || die "Could not open $outfn\n";

	# Read/Process fasta file, extracting records as we need them
	my $offset=0;
	print STDERR "Processing FASTA file...\n";
	my ($defline, $prev_defline, $sequence);
	while(<$infh>){
		chomp;
		
		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, $sequence, $offset, $select_arr, $outfh);
				$sequence="";
				$offset++;
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
		
	}
	process_record($prev_defline, $sequence, $offset, $select_arr, $outfh);

	close($infh);
	close($outfh);

	print STDERR "Done.\n\n";
}

print STDERR "Completed.\n";

###############################################################################
###############################################################################

sub select_samples{
	my $num_samples=shift;
	my $num_records=shift;

	my @selected_arr;

	# If num samples requested is less than half
	# of what is available, randomly pick which samples to keep	
	#
	# If num samples needed is more than half of 
	# of the available, randomly pick what to remove

	if(1.0*$num_samples/$num_records < .5){

		for(my $i=0; $i<$num_records; $i++){
			$selected_arr[$i]=0;
		}

		print STDERR "Picking to include...\n";
		for(my $p=0; $p<$num_samples;){
			my $offset=int(rand($num_records));
			if($selected_arr[$offset]==0){
				$selected_arr[$offset]=1;
				$p++;
			}
		}
	}else{

		for(my $i=0; $i<$num_records; $i++){
			$selected_arr[$i]=1;
		}

		print STDERR "Picking to exclude...\n";
		for(my $p=$num_records; $p>$num_samples;){
			my $offset=int(rand($num_records));
			if($selected_arr[$offset]==1){
				$selected_arr[$offset]=0;
				$p--;
			}
		}
	}

	return(\@selected_arr);
}

sub count_records{
	my $name=shift;
	open(FH, "<$name") || die "Could not open $name\n";
	
	my $num_records=0;
	while(<FH>){
		if($_=~/^>/){
			$num_records++;
		}
	}
	
	close(FH);
	return($num_records);
}

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;
	my $offset = shift;
	my $select_arr= shift;
	my $fh = shift;

	if($defline eq "" && $sequence eq ""){
		print STDERR "Empty Record.\n";
		return;
	}

	if(!defined($select_arr) || ${$select_arr}[$offset]){
		# Output record if the selection array is defined
		# or it the record has been selected

		print $fh "$defline\n";

		my $length=length($sequence);
		my $width=80;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print $fh substr($sequence, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);
	}
}

#------------------------------------------------------------------------------

