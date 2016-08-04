#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_s $opt_o);
getopts("i:s:o:");


my $usage = "usage:
$0
	-i <input summary table>
	-s <sample list to extract>
	-o <output summary table root>

This script will read in the sample list and extract the
samples that are specified in the list.  It will exclude categories
that will be 0 in the extracted list.

To save memory, the script will read through the
input summary table twice.  The first time to determine
which categories are non-zero in the extracted summary
table, and then the second time to write out the new
summary table while excluding categories that are all
zeros in the sample subset.

";

if(!(
	defined($opt_i) &&
	defined($opt_s) &&
	defined($opt_o)
)){
	die $usage;
}


my $input_st=$opt_i;
my @ext=(".summary_table.tsv", ".summary_table.xls");
my ($name, $path, $extension)=fileparse($input_st, @ext);

my $samp_list=$opt_s;
my $output_fn=$opt_o;

$output_fn=~s/\.summary_table\.tsv//;
$output_fn="$output_fn\.summary_table.tsv";

print STDERR "\n";
print STDERR "Input Summary Table: $input_st\n";
print STDERR "Sample Keep List: $samp_list\n";
print STDERR "Output Summary Table: $output_fn\n";
print STDERR "\n";

###############################################################################

sub load_list{
	my $fn=shift;
	my %list_hash;
	
	open(FH, "<$fn") || die "Could not open $fn\n";
	while(<FH>){
		chomp;
		my ($item, $junk)=split "\t", $_;
		$list_hash{$item}=1;
	}
	close(FH);
	
	return(\%list_hash);
}

###############################################################################

print STDERR "Loading Sample List...\n";
my $list_hash_ref=load_list($samp_list);

###############################################################################
# First pass, identify which categories are non-zero

my %original_totals;
print STDERR "Skimming summary table for targeted samples...\n";
open(FH, "<$input_st") || die "Could not open $input_st\n";

# Read in header info
my $hdr_line=<FH>;
chomp $hdr_line;

my @hdr_col=split "\t", $hdr_line;

my @categories=@hdr_col;
shift @categories;	# Shift out sample_id
shift @categories;	# Shift out total

my $num_categories=$#categories+1;
print STDERR "Num categories found: $num_categories\n";

# Initialize counts
my @cumulative_counts;
for(my $i=0; $i<$num_categories; $i++){
	$cumulative_counts[$i]=0;
}

# Read in counts, but process 1 line at a time.
while(<FH>){
	chomp;

	my @data_col=split "\t", $_;

	my $sample_id=shift @data_col;
	my $total=shift @data_col;

	# Count up categories across kept samples
	if(defined(${$list_hash_ref}{$sample_id})){
		print STDERR "$sample_id found.\n";
		my $samp_tot=0;
		for(my $i=0; $i<$num_categories; $i++){
			$cumulative_counts[$i]+=$data_col[$i];
			$samp_tot+=$data_col[$i];
		}
		$original_totals{$sample_id}=$samp_tot;
	}
}

close(FH);

###############################################################################

# Identify indices of non-zero categories
my @non_zero;
for(my $i=0; $i<$num_categories; $i++){
	if($cumulative_counts[$i]>0){
		push @non_zero, $i;
	}
}

my $num_non_zero=$#non_zero+1;
print STDERR "Number of non-zero categories found in subset: $num_non_zero\n";

###############################################################################
# Second pass, write out samples and categories that were targeted

open(FH, "<$input_st") || die "Could not open $input_st\n";
open(FH_OUT, ">$output_fn") || die "Could not open $output_fn\n";

# Read in header info
my $hdr_line=<FH>;
chomp $hdr_line;
my @hdr_col=split "\t", $hdr_line;

my @categories=@hdr_col;
my $samp_id_str=shift @categories;      # Shift out sample_id
my $total_str=shift @categories;      # Shift out total

# Write out categories that will be kept
print FH_OUT "$samp_id_str\t$total_str\t";
my @new_cat;
for(my $i=0; $i<$num_non_zero; $i++){
	push @new_cat, $categories[$non_zero[$i]];	
}
print FH_OUT (join "\t", @new_cat) . "\n";

# Write out samples/counts for kept categories
while(<FH>){
        chomp;

        my @data_col=split "\t", $_;

        my $sample_id=shift @data_col;
        my $total=shift @data_col;

        if(defined(${$list_hash_ref}{$sample_id})){
		print FH_OUT "$sample_id\t$total\t";
		my @counts;
		# Write out counts that are non-zero across all subset samples
		my $samp_tot=0;
                for(my $i=0; $i<$num_non_zero; $i++){
			my $ct=$data_col[$non_zero[$i]];	
			push @counts, $ct;
			$samp_tot+=$ct;
                }
		if($samp_tot!=$original_totals{$sample_id}){
			print STDERR "ERROR: Written totals ($samp_tot) not equal to original ($original_totals{$sample_id})\n";
		}
		print FH_OUT (join "\t", @counts) . "\n";
        }
}

close(FH);
close(FH_OUT);

###############################################################################

print STDERR "Done.\n";




