#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_l $opt_M $opt_w $opt_b $opt_e);

my $WIDTH=60;
my $DEFAULT_BEGIN=1;
my $DEFAULT_END=2;

getopts("f:l:M:w:b:e:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-l <list of trim coordinates>
	[-M <mask char/quality value>]
	[-w <width, default $WIDTH>]

	Use alternate trim columns, counting from 0.
	[-b <clear begin column, default=$DEFAULT_BEGIN>]
	[-e <clear end column, default=$DEFAULT_END>]

	Given a list of clear range:

	<seq_id>\\t<begin>\\t<end>\\n

	Will trim the input fasta file.  Input maybe sequence or quality.
	Deflines are preserved.  If -M flag is used, then sequences are masked
	with the specified character.  If the input is a qual file, the
	mask will be a quality value.

	By default, the second and third column will be the begin and end
	used, respectively.  But you can specify a different -b/-e if the input
	file has a different format.

	Output goes to STDOUT.

";

if(!(
	defined($opt_f) && 
	defined($opt_l))){
	die $usage;
}

my $mask_char=undef;
if(defined($opt_M)){
	$mask_char=$opt_M;
}

my $width=$WIDTH;
if(defined($opt_w)){
	$width=$opt_w;
}

my $begin_col=$DEFAULT_BEGIN;
if(defined($opt_b)){
	$begin_col=$opt_b;
}

my $end_col=$DEFAULT_END;
if(defined($opt_e)){
	$end_col=$opt_e;
}


###############################################################################
# Read in trim list

print STDERR "Loading trim list: $opt_l\n";
open(LIST_FH, "<$opt_l") || die "Could not open $opt_l\n";

my %list_hash;
my $list_length=0;
while(<LIST_FH>){
	chomp;

	my @in=split /\t/, $_;
	my $id=$in[0];
	my $begin=$in[$begin_col];
	my $end=$in[$end_col];

	$list_hash{$id}="$begin#$end\n";
	$list_length++;
}

close(LIST_FH);
print STDERR "Done. $list_length regions loaded.\n";

###############################################################################
# Read in features

my $num_found=0;

my $UNKNOWN=-1;
my $QUAL=2;
my $RESIDUE=1;

my $seq_type=$UNKNOWN;

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

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
		if($seq_type==$UNKNOWN){
			my $sample=substr($_,0,10);
			$sample=~s/\s+//g;
			if($sample=~/\d/){
				$seq_type=$QUAL;
				print STDERR "FASTA type is Quality.\n";	
			}elsif($sample=~/[A-Za-z]/){
				$seq_type=$RESIDUE;
				print STDERR "FASTA type is Sequence.\n";	
			}else{
				die "Could not determine sequence type.\n";
			}
		}

		if($seq_type==$QUAL){
			$sequence.=" " . $_;
		}else{
			$sequence.=$_;
		}
	}
}
process_record($prev_defline, $sequence);

close(FASTA_FH);

print STDERR "$num_found out of $list_length found\n";
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse id from defline: $defline\n";
	}

	$sequence=~s/^\s+//;
	$sequence=~s/\s+$//;

	my ($begin, $end)=(-1,-1);
	if(defined($list_hash{$id})){
		($begin, $end)=split /#/, $list_hash{$id};
		$num_found++;
        # skip the 0, 0 ones
        return if $begin == 0 && $end == 0;
	}
	
	print STDOUT "$defline\n";

	if($seq_type==$QUAL){
		my @qual=split /\s+/, $sequence;
		dump_qual(\@qual, $begin, $end);
	}elsif($seq_type==$RESIDUE){
		$sequence=~s/\s+//g;
		dump_sequence(\$sequence, $begin, $end);
	}
}

#------------------------------------------------------------------------------
	
sub dump_sequence{
	my $seq_ref=shift;
	my $begin=shift;
	my $end=shift;

	my $length=length(${$seq_ref});
	
	if($begin!=-1){
		if(defined($mask_char)){
			#mask
			substr(${$seq_ref}, $end, ($length-$end))=($mask_char x ($length-$end));
			substr(${$seq_ref}, 0, $begin)=($mask_char x $begin);
		}else{
			#trim
			${$seq_ref}=substr(${$seq_ref}, $begin, $end-$begin);
		}

		$length=length(${$seq_ref});
	}

	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr(${$seq_ref}, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
		
}

#------------------------------------------------------------------------------

sub dump_qual{
	my $qual_arr_ref=shift;
	my $begin=shift;
	my $end=shift;

	my $length=$#{$qual_arr_ref}+1;
	
	my @arr;
	if($begin!=-1){
		@arr=splice(@{$qual_arr_ref}, $begin, ($end-$begin));
	}else{
		@arr=@{$qual_arr_ref};
	}

	if(defined($mask_char)){
		my (@front, @end);
		for(my $i=0; $i<$begin; $i++){
			push @front, $mask_char;
		}
		for(my $i=$end; $i<$length; $i++){
			push @end, $mask_char;
		}
		unshift @arr, @front;
		push @arr, @end;
	}
	
	my @line;
	my $remaining=$#arr+1;

	while($remaining>0){
		@line=splice(@arr, 0, $width);
		my $outline=join " ", @line;
		print STDOUT "$outline\n";
		$remaining=$#arr+1;
	}

}

