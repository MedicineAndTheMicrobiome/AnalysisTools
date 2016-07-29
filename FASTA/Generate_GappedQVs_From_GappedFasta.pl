#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_g $opt_q $opt_c $opt_o $opt_a);

my $LN10=log(10);

getopts("g:q:c:o:a");
my $usage = "usage: 
$0 
	-g <gapped fasta file>
	-q <ungapped quality value file>
	[-c <read used range>]
	-o <gapped quality file, output>

	[-a (flag to fill gaps with average of QV around it)]

	This program will read in a gapped fasta file, and use it to
	introduce gaps into a QV file.  If the -c option is not used
	then the number of quality values better be equal to the number of
	nongaps in the fasta file, AND be in the same orientation.

	if the -a flag is used, if there is a gap, then the average of the QV of
	the bases surrounding it will be used.  Gaps in the beginning/ending of
	the sequence will remain as -'s.  


";

if(
	!(defined($opt_g)) ||
	!(defined($opt_q)) ||
	!(defined($opt_o))
){
	die $usage;
}

my $gapped_fasta=$opt_g;
my $quality_fasta=$opt_q;
my $read_ranges_file=$opt_c;
my $gapped_qual=$opt_o;
my $use_avg_for_gaps=defined($opt_a);

my %quality_hash;
my %sequence_hash;
my @id_arr;
my %ranges_hash;

open(FASTA_FH, "<$gapped_fasta") || die "Could not open $gapped_fasta\n";
open(QUAL_FH, "<$quality_fasta") || die "Could not open $quality_fasta\n";
open(GAP_QUAL_FH, ">$gapped_qual") || die "Could not open $gapped_qual, for writing\n";

###############################################################################

if(defined($read_ranges_file)){
	open(RANGES_FH, "<$read_ranges_file") || die "Could not open $read_ranges_file\n";
	while(<RANGES_FH>){
		chomp;
		my ($id, $begin, $end, $ori)=split /\t/, $_;
		$ranges_hash{$id}="$begin#$end#$ori";
	}
	close(RANGES_FH);
}

###############################################################################

print STDERR "Processing Gapped FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_sequence($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_sequence($prev_defline, $sequence);

#-----------------------------------------------------------------------------

print STDERR "Processing Quality FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<QUAL_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_quality($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=($_ . " ");
	}
}
process_quality($prev_defline, $sequence);

###############################################################################

my $num_gapped_sequences=$#id_arr+1;
print STDERR "Num gapped sequences: $num_gapped_sequences\n";

foreach my $id(@id_arr){
	#print "$id\n";
	#print $quality_hash{$id}\n";
	#print "$sequence_hash{$id}\n\n";
	
	if(!defined($quality_hash{$id})){
		die "Could not find quality values for '$id'\n";
	}

	my ($begin, $end, $ori);
	if(defined($ranges_hash{$id})){
		($begin, $end, $ori)=split /#/, $ranges_hash{$id};
	}else{
		($begin, $end, $ori)=(0, $#{$quality_hash{$id}}+1, 1);
	}
	my $gapped_quality_arr_ref=gap_quality(\$sequence_hash{$id}, $quality_hash{$id}, $begin, $end, $ori);

	if($use_avg_for_gaps){
		$gapped_quality_arr_ref=use_avg_for_gap_quality($gapped_quality_arr_ref);
	}

	output_gapped_quality(\*GAP_QUAL_FH, $id, $gapped_quality_arr_ref, \$sequence_hash{$id});
}

print STDERR "Completed.\n";


###############################################################################

sub use_avg_for_gap_quality{
	my $inp_qual_ref=shift;

	my @new_gap_qual=@{$inp_qual_ref};

	my $len=$#{$inp_qual_ref}+1;

	my $left_pad=0;
	my $right_pad=$len-1;

	# Compute where gaps are starting
	while(${$inp_qual_ref}[$left_pad] eq "-"){
		$left_pad++;
	}

	# Compute where gaps are ending
	while(${$inp_qual_ref}[$right_pad] eq "-"){
		$right_pad--;
	}

	# Go through gaps, and compute average QV based on bases around them
	for(my $i=$left_pad+1; $i<$right_pad; $i++){
		if(${$inp_qual_ref}[$i] eq "-"){
			my ($l, $r)=find_non_gap_around_pos($i, $inp_qual_ref);
			$new_gap_qual[$i]=sprintf("%3.2f", avg_qv(${$inp_qual_ref}[$l], ${$inp_qual_ref}[$r]));
		}
	}

	return(\@new_gap_qual);

}

#-------------------------------------------------------------------------------
# Code to make sure we average QV properly

sub qv_to_errPerBase{
	my $qv=shift;
	return(exp(-$LN10*$qv/10));
}

sub errPerBase_to_qv{
	my $errs=shift;
	return(-10*log($errs)/$LN10);
}

sub avg_qv{
	my $a=shift;
	my $b=shift;
	my $epb_a=qv_to_errPerBase($a);
	my $epb_b=qv_to_errPerBase($b);
	return(errPerBase_to_qv(($epb_a+$epb_b)/2.0));
}

#-------------------------------------------------------------------------------

sub find_non_gap_around_pos{
	my $pos=shift;
	my $qual_arr_ref=shift;

	# Given a position and the gapped quality array, it will return the first non gap positions
	# up and downstream of the position specified.

	#look up
	my $up=$pos-1;
	while(${$qual_arr_ref}[$up] eq "-"){
		$up--;
	}

	#look down
	my $down=$pos+1;
	while(${$qual_arr_ref}[$down] eq "-"){
		$down++;
	}
	
	return($up, $down);

}

################################################################################

sub gap_quality{
	my $seq_ref=shift;
	my $qual_arr_ref=shift;
	my $begin=shift;
	my $end=shift;
	my $ori=shift;

	my $ungapped_seq=${$seq_ref};
	$ungapped_seq=~s/-//g;
	my $ungapped_seq_len=length($ungapped_seq);
	my $used_qual_len=($end-$begin);

	#print "  ungapped seq length = $ungapped_seq_len\n";
	#print "  used quality length = $used_qual_len\n";
	if($ungapped_seq_len != $used_qual_len){
		die "Error: mismatch between ungapped sequence length and used quality value length: $ungapped_seq_len != $used_qual_len\n";
	}

	my @nuc_arr=split //, ${$seq_ref};	
	my @gapped_qual;

	my $q_i;
	my $inc;

	# Determine which direction to start aligning quality from.
	if($ori==1){
		$ori=$begin;
		$inc=1;
	}else{
		$ori=$end-1;
		$inc=-1;
	}

	for(my $i=0; $i<=$#nuc_arr; $i++){
		#print "$nuc_arr[$i]\n";
		if($nuc_arr[$i] ne "-"){
			push @gapped_qual, ${$qual_arr_ref}[$q_i];
			$q_i+=$inc;
		}else{
			push @gapped_qual, "-";
		}
	}

	return(\@gapped_qual);
}

#-------------------------------------------------------------------------------

sub output_gapped_quality{
	my $fh=shift;
	my $id=shift;
	my $gapped_qual_arr_ref=shift;
	#my $sequence_arr_ref=shift;

	print {$fh} ">$id\n";

	my $arr_len=$#{$gapped_qual_arr_ref}+1;
	my $width=50;	
	my $offset;

        my $abs_pos=0;
        while($abs_pos<$arr_len){
                for(my $i=0; $i<$width && $abs_pos<$arr_len; $i++){
                        if($i>0){
                                print {$fh} " ${$gapped_qual_arr_ref}[$abs_pos]";
                        }else{
                                print {$fh} "${$gapped_qual_arr_ref}[$abs_pos]";
                        }
                        $abs_pos++;
                }
                print {$fh} "\n";
        }
}

###############################################################################
###############################################################################

sub process_quality{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;	
	}

	$sequence=~s/^\s+//;
	$sequence=~s/\s+$//;
	my @qv=split /\s+/, $sequence;

	$quality_hash{$id}=\@qv;

}

#------------------------------------------------------------------------------

sub process_sequence{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;	
	}

	push @id_arr, $id;
	$sequence_hash{$id}=$sequence;
}

###############################################################################




