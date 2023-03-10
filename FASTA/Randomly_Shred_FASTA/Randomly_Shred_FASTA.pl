#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_s $opt_n $opt_r $opt_R);

getopts("f:s:n:r:R");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-n <Num fragments>
	-m <mean length>
	-s <stdev length>
	[-c circularize genome]

	This program will take the input fasta file and randomly generate
	a set of reads with the specified mean and standard deviation.
	If you want the genome to be circularized, use the -c flag or else
	you may end up with many short sequences off the end of the sequence.

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

if($num_records<=$sample_size && !$with_replacement){
	die "Num records is less than sample size. You can't sample without replacement.\n";
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

