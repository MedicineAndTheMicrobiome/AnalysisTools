#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use File::Basename;
use Cwd;
use vars qw($opt_f $opt_t $opt_c $opt_n $opt_m $opt_o $opt_d $opt_w);

my $WIDTH=50;
my $FASTA_TYPE="seq";

getopts("f:t:c:m:n:o:d:w:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-t <Tab Separated Value File to Split FASTA records by>
	-c <columns to separate by, as a list of numbers separated by commas, first column is 0>
	[-n <FASTA file type, 'seq' or 'qual', default='seq'>]
	[-m <id mapping file, for example barcode-to-library name>]
	[-w <width, default=$WIDTH>]
	[-o <Output Root FASTA Filename>]
	[-d <Output Directory>]

	Reads in the Tab Separated Value File, and based on the columns specified,
	will separate out a set a fasta files based on those keys.

	The TSV file must have the first column (0), be the id that is used in the
	fasta file.

	If the FASTA file contains the same ID twice, it is eliminated the second
	time it is encountered.

	The id mapping file can be specified to rename the identifiers in the specified column
	to a new name.

";

if(!(
	defined($opt_f) && 
	defined($opt_t) && 
	defined($opt_c))){
	die $usage;
}

my $in_fasta=$opt_f;
my $in_tsv=$opt_t;
my @columns=split ",", $opt_c;
my $output_dir;
my $width;
my $mapping_file=$opt_m;
my $out_fasta_root;
my $id_map_hash_ref;
my $fasta_type=$FASTA_TYPE;

###############################################################################

# Width
if(!defined($opt_w)){
	$width=$WIDTH;
}else{
	$width=$opt_w;
}

# Output Dir
if(defined($opt_d)){
	$output_dir=$opt_d;
	if(!(-e $output_dir)){
		`mkdir $output_dir`;	
	}
}else{
	$output_dir=cwd();
}

# Output root filename
if(!defined($opt_o)){
	$out_fasta_root=$in_fasta;
}else{
	$out_fasta_root=$opt_o;
}

# Output path
if(defined($opt_d)){
	($out_fasta_root)=fileparse($out_fasta_root);
	$out_fasta_root="$output_dir/$out_fasta_root"
}

# FASTA type
if(defined($opt_n)){
	$fasta_type=$opt_n;
}

if($fasta_type ne "seq" && $fasta_type ne "qual"){
	die "Invalid FASTA type: $fasta_type\n";
}

# Mapping file
if(!defined($mapping_file)){
	$mapping_file="";
}

print STDERR "\n";
print STDERR "Column Number(s) to separate sequences by: " . (join "/", @columns) ."\n";
print STDERR "Output Filename root: $out_fasta_root\n";
print STDERR "FASTA Type: $fasta_type\n";
print STDERR "Residues/line: $width\n";

if($mapping_file ne ""){
	print STDERR "Mapping file: $mapping_file\n";
	$id_map_hash_ref=load_id_mapping($mapping_file);
}
print STDERR "\n";

###############################################################################

sub load_id_mapping{
	my $map_fname=shift;
	my %id_hash;
	print STDERR "Loading ID mapping file: $map_fname\n";
	open(FH, "<$map_fname") || die "Could not open $map_fname\n";
	while(<FH>){
		chomp;
		my ($old, $new)=split "\t", $_;
		$id_hash{$old}=$new;
	}
	close(FH);
	my $num_uniq_ids=keys %id_hash;
	print STDERR "   Num IDs in map: $num_uniq_ids\n";
	return(\%id_hash);
}

###############################################################################
# Build hash of where each sequence should be partitioned into

my %seq_mapping_hash;
my %partitions_hash;
open(TSV_FH, "<$in_tsv") || die "Could not open $in_tsv for reading.\n";

print STDERR "Loading $in_tsv...\n";
while(<TSV_FH>){
	chomp;
	my @cols=split /\t/, $_;
	
	# Grab the columns we need
	my @key_comp;
	foreach my $targ_col(@columns){
		push @key_comp, @cols[$targ_col];
	}

	# Join the columns we grabbed into a key
	my $key=join ".", @key_comp;
	$seq_mapping_hash{$cols[0]}=$key;
	$partitions_hash{$key}=1;
}
print STDERR "ok.\n\n";

close(TSV_FH);

my $num_seqs_described=keys %seq_mapping_hash;
my $num_partitions=keys %partitions_hash;
print STDERR "Split Table Information:\n";
print STDERR "   Num sequences described: $num_seqs_described\n";
print STDERR "   Num unique partitions:   $num_partitions\n";
print STDERR "\n";

###############################################################################
# Make sure files open before wasting any time doing anything

my $extension;
if($fasta_type eq "seq"){
	$extension="fasta";
}else{
	$extension="qual";
}

foreach my $partition(sort keys %partitions_hash){
	
	my $partition_name=$partition;
	if($mapping_file ne ""){
		$partition_name=${$id_map_hash_ref}{$partition};
		if(!defined($partition_name)){
			$partition_name=$partition;
		}
	}

	$partitions_hash{$partition}=new FileHandle ">$out_fasta_root\.$partition_name\.$extension";
}

###############################################################################
# Read in features

my $unmapped=0;

open(FASTA_FH, "<$in_fasta") || die "Could not open $in_fasta\n";
print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline);
my @sequence;
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($#sequence!=-1){
			if($fasta_type eq "seq"){
				my $test_char=$sequence[$#sequence];
				if($test_char=~/[0-9]/){
					print STDERR "Error:  Numbers found in sequence FASTA file.  Are you sure this is not a quality file?\n";
					print STDERR "$prev_defline\n";
					print STDERR (join " ", @sequence) . "\n";
					die "Invalid characters in FASTA sequence file.\n";
				}
			}
			process_record($prev_defline, \@sequence);
			@sequence=();
		}
		$prev_defline=$defline;
	}else{
		if($fasta_type eq "seq"){
			push @sequence, (split //, $_);
		}else{
			push @sequence, (split /\s+/, $_);
		}
	}
}
process_record($prev_defline, \@sequence);

close(FASTA_FH);

print STDERR "Num Sequences without a Home: $unmapped\n";
print STDERR "Completed.\n";

###############################################################################

my %outputted_seq_hash;

sub process_record{
	my $defline = shift;
	my $sequence_ref = shift;
	
	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}

	if(defined($seq_mapping_hash{$id})){
	
		if(!defined($outputted_seq_hash{$id})){
		
			if(!defined($seq_mapping_hash{$id}) || !defined($partitions_hash{$seq_mapping_hash{$id}})){
				print STDERR "Error: (are you sure you picked the right columns for input?)\n";
				print STDERR "  id=$id\n";
				print STDERR "  seq_mapping_hash=$seq_mapping_hash{$id}\n";
				print STDERR "  partitions_hash=$partitions_hash{$seq_mapping_hash{$id}}\n";
			}
			print {$partitions_hash{$seq_mapping_hash{$id}}} "$defline\n";

			my $length=$#sequence+1;
			my $pos=0;
			do{
				my $out_width=($width>$length)?$length:$width;

				my @out_arr=splice @{$sequence_ref}, 0, $out_width;
	
				if($fasta_type eq "seq"){
					print {$partitions_hash{$seq_mapping_hash{$id}}} (join "", @out_arr) . "\n";
				}else{
					print {$partitions_hash{$seq_mapping_hash{$id}}} (join " ", @out_arr) . "\n";
				}

				$length-=$width;
			}while($length>0);
			$outputted_seq_hash{$id}=1;
		}

	}else{
		$unmapped++;
	}

}

#------------------------------------------------------------------------------

