#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_q $opt_o);

my $OFFSET=33;

getopts("f:q:o:");
my $usage = "
usage: 
$0 
	-f <input fasta sequence file>
	-q <input fasta quality file>
	[-o <output fastq file>]

	This script will read in both sequence and quality file
	and generate a fastq file.

	It will confirm that the number of residues match the
	number of quality values and that the id's between the
	two files match up, and the quality values are within range.

	You should confirm that the sequences are in the same
	order, or else an error will be thrown.

	The encoding  offset will be +$OFFSET, and the expected range will 
	be between 0 and 41.

";

if(!(
	defined($opt_f) &&
	defined($opt_q))){
	die $usage;
}

my $fasta_name=$opt_f;
my $qual_name=$opt_q;
my $fastq_name=$opt_o;

if(!defined($fastq_name)){
	$fastq_name=$fasta_name;
	$fastq_name=~s/\.fasta$//;
	$fastq_name.=".fastq";
}

###############################################################################

print STDERR "Input Sequence FASTA: $fasta_name\n";
print STDERR "Input Quality FASTA: $qual_name\n";
print STDERR "Output FASTQ: $fastq_name\n";


open(SEQ_FH,  "<$fasta_name") || die "Could not open $fasta_name for reading.\n";
open(QUAL_FH, "<$qual_name") || die "Could not open $qual_name for reading.\n";
open(FASTQ_FH, ">$fastq_name") || die "Could not open $fastq_name for writing.\n";

###############################################################################

# Static Variables
my $qual_next_id;
my $seq_next_id;

sub read_fasta_qual_record{
	my $fh=shift;

	my $qual_buffer;
	while(<$fh>){
		chomp;
		my $line=$_;
		
		if($line=~/^>(\S+)/){
			my $cur_id=$qual_next_id;
			$qual_next_id=$1;

			if($qual_buffer ne ""){
				$qual_buffer=~s/^\s+//;
				$qual_buffer=~s/\s+$//;
				my @qual_val=split /\s+/, $qual_buffer;
				return($cur_id, \@qual_val);
			}
		}else{
			$qual_buffer.= (" " . $line);	
		}
	}

	if(eof($fh)){
		$qual_buffer=~s/^\s+//;
		$qual_buffer=~s/\s+$//;
		my @qual_val=split /\s+/, $qual_buffer;
		return($qual_next_id, \@qual_val);
	}
	
}

#------------------------------------------------------------------------------

sub read_fasta_sequence_record{
	my $fh=shift;

	my $seq_buffer;
	while(<$fh>){
		chomp;
		my $line=$_;
		
		if($line=~/^>(\S+)/){
			my $cur_id=$seq_next_id;
			$seq_next_id=$1;

			if($seq_buffer ne ""){
				return($cur_id, $seq_buffer);
			}
		}else{
			$seq_buffer.= $line;	
		}
	}

	if(eof($fh)){
		return($seq_next_id, $seq_buffer);
	}
	
}

#------------------------------------------------------------------------------

sub write_fastq{
	my $id=shift;
	my $sequence=shift;
	my $qual_arr_ref=shift;
	my $fh=shift;

	# Confirm lengths are the same
	my $seq_length=length($sequence);
	my $qual_length=$#{$qual_arr_ref}+1;

	if($seq_length != $qual_length){
		die "Error: Sequence and Quality lengths don't match: $seq_length (seq) != $qual_length (qual).\n";
	}
	
	# Encode in string
	my @qual_char;
	foreach my $qual_val(@{$qual_arr_ref}){
		if($qual_val < 0 || $qual_val > 41){
			die "Error: Quality value ($qual_val), out of range (0-41)\n";
		}
		push @qual_char, chr($qual_val+$OFFSET);
	}
	my $qual_code=join "", @qual_char;

	# Output Fastq record
	print {$fh} "\@$id\n";
	print {$fh} "$sequence\n";
	print {$fh} "+\n";
	print {$fh} "$qual_code\n";

	return;
}

###############################################################################

print STDERR "\nReading in \n\t$fasta_name\nand\n\t$qual_name\nrecords in parallel...\n\n";

my $endoffile=0;
my $num_records_proc=0;
while(!$endoffile){

	my ($qual_id, $qual_val_arr_ref)=read_fasta_qual_record(\*QUAL_FH);
	my ($seq_id, $sequence)=read_fasta_sequence_record(\*SEQ_FH);

	#print STDERR "'$seq_id' / '$qual_id'\n";
	#print STDERR "$sequence\n";
	#print STDERR (join " ", @{$qual_val_arr_ref}) . "\n";
	#print STDERR "\n";

	if($seq_id ne $qual_id){
		die "Error!  $seq_id ne $qual_id.  Please confirm records are in the same order.\n";
	}

	write_fastq($seq_id, $sequence, $qual_val_arr_ref, \*FASTQ_FH);
	$num_records_proc++;
	if(!($num_records_proc%10000)){
		print STDERR ".";
	}

	$endoffile=(eof(SEQ_FH) || eof(QUAL_FH));
}
print STDERR "\n";

if(!eof(SEQ_FH) || !eof(QUAL_FH)){
	die "Error the end of file not reached for both sequence and quality file!!!\n";
}

print STDERR "Num Records Processed: $num_records_proc\n";
print STDERR "Completed.\n";

###############################################################################
