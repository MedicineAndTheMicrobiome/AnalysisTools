#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_s $opt_n $opt_r $opt_R $opt_L $opt_S);

getopts("f:s:n:r:RLS:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-n <sample size, n>
	-s <number of samples (bootstraps)>
	-r <output root>
	[-R (with replacement flag)]
	[-L (Limit sample size to number of sequeces in FASTA)]
	[-S <Set random number seed, default is not set so PERL uses time>]

	This program will randomly sample s sequences from the fasta file f,
	n number of times.  Output will go to:

	<output root>.<index>.fasta

	Sampling will be with replacement and uniformly distributed.

";

if(!(
	defined($opt_f) && 
	defined($opt_s) && 
	defined($opt_n) && 
	defined($opt_r))){
	die $usage;
}

###############################################################################

my $fasta=$opt_f;
my $num_samples=$opt_s;
my $sample_size=$opt_n;
my $output_root=$opt_r;
my $with_replacement=defined($opt_R);
my $limit_by_num_reads=defined($opt_L);
my $random_number_seed=undef;

if(defined($opt_S)){
	$random_number_seed=$opt_S;
	srand($random_number_seed);
	print STDERR "Random number seed: $random_number_seed\n";
}

print STDERR "Input fasta: $fasta\n";
print STDERR "Number of samples: $num_samples\n";
print STDERR "Sample size: $sample_size\n";
print STDERR "Output root: $output_root\n";

# Open output fasta files
my @fh;
my $digits=int(log($num_samples)/log(10))+1;
for(my $i=0; $i<$num_samples; $i++){
	$fh[$i]=FileHandle->new();
	my $offset=sprintf("%0$digits" . "i", $i);
	my $fname="$output_root\.$offset\.fasta";
	$fh[$i]->open(">$fname") || die "Could not open $fname\n";
}

# Determine population size
print STDERR "Counting number of records in FASTA file\n";
my $num_records=count_records($fasta);
print STDERR "Number of records: $num_records\n";

if($num_records<$sample_size){
	print STDERR "Requested sample size less than number of available sequences\n";

	if($limit_by_num_reads){
		print STDERR "Taking option to limit number of sequences in sample.\n";
		$sample_size=$num_records;
	}elsif(!$with_replacement){
		print STDERR "Num records is less than sample size. You can't sample without replacement.\n";
		print STDERR "Either limit the number of sequences per sample, or allow sampling with replacent.\n";
		die;
	}else{
		print STDERR "Sampling with replacement.\n";
	}
}

# Determine which records should go in which file, randomly
my $subset_hash_ref=compute_samples($num_samples, $sample_size, $num_records, $with_replacement);

# Read/Process fasta file, extracting records as we need them
open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

my $offset=0;
print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
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

close(FASTA_FH);

print STDERR "Completed.\n";

###############################################################################
###############################################################################

sub compute_samples{
	my $num_samples=shift;
	my $sample_size=shift;
	my $population_size=shift;
	my $with_replacement=shift;

	my %sample_hash;

	for(my $s=0; $s<$num_samples; $s++){

		my %selected_hash;

		for(my $p=0; $p<$sample_size;){
			my $offset=int(rand($population_size));
			#print "$offset: $s\n";
			if($with_replacement || (!$with_replacement && !defined($selected_hash{$offset}))){
				push @{$sample_hash{$offset}}, $s;
				$selected_hash{$offset}=1;
				$p++;
			}
		}
	}
	return(\%sample_hash);
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

	if(defined(${$subset_hash_ref}{$offset})){
		
		my $samples_ref=${$subset_hash_ref}{$offset};

		foreach my $sample(@{$samples_ref}){

			print {$fh[$sample]} "$defline\n";

			my $length=length($sequence);
			my $width=70;
			my $pos=0;
			do{
				my $out_width=($width>$length)?$length:$width;
				print {$fh[$sample]} substr($sequence, $pos, $width) . "\n";
				$pos+=$width;
				$length-=$width;
			}while($length>0);
		}
	}

	$offset++;
}

#------------------------------------------------------------------------------

