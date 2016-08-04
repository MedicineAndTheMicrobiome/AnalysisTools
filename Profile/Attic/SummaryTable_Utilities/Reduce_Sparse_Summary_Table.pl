#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_o $opt_t $opt_p $opt_m);
getopts("i:o:t:p:m:");


my $usage = "usage:
$0
	-i <input summary table>
	[-o <output root>]

	Keep criteria:
	[-t <number of top categories to keep, e.g. 200>]
	[-p <percentage of top categories to keep, e.g. 95>]
	[-m <minimum counts to keep, e.g. 2>]

	This script will read in a summary table that is very sparse (i.e.
	mostly 0's) and then try to reduce it by one of three methods.
	However, first, each sample is normalized to account for differences in
	sequencing depth.  Then across all samples, the average abundance
	for each category is computed.  Categories are then sorted
	by average abundance.  (If you have negative controls where
	the categories are low and spurious, you may want to remove
	these for a final analysis, because they may skew the average
	abundances across all your samples.)

	The keep criteria are:

	Top N categories:
	If the -t option is taken, the top N most abundant catgories
	are saved, the remaining are placed into a Remaining category.
	
	Top P Percentage:
	If the -p option is taken, then top most abundant categories
	that sum up to the P% of the abundances are saved, and the
	remaining are placed into Remaining category.

	Minimum M Counts:
	If the -m option is taken, then a category is removed if
	across all samples the count is less than M and there are
	no cases for which a category with count less than M has
	an abundance exceeding it.  This is a little confusing.  
	Let's say you want to remove singletons.  A singleton
	may not be low abundance if the total number of reads are
	low for that sample.  We only want to remove singletons
	that are associated with high sequencing depth.  Using this
	convoluted criteria, if a doubleton exists at an abundance
	A, then no singletons with an abundance >= A will be removed.

";

if(!(
	defined($opt_i)
)){
	die $usage;
}

###############################################################################

my $input_st=$opt_i;
my @ext=(".summary_table.tsv", ".summary_table.xls");
my ($name, $path, $extension)=fileparse($input_st, @ext);

my $output_root="$path/$name";

if(defined($opt_o)){
	$output_root=$opt_o;
}

my $keep_top_N=0;
my $keep_top_P=0;
my $keep_min_M=0;
if(defined($opt_t)){
	$keep_top_N=$opt_t;
}
if(defined($opt_p)){
	$keep_top_P=$opt_p;
}
if(defined($opt_m)){
	$keep_min_M=$opt_m;
}

print STDERR "Input File: $input_st\n";
print STDERR "Output Root: $output_root\n";

if($keep_top_N){
	print STDERR "Keeping top $keep_top_N categories.\n";
}
elsif($keep_top_P){
	print STDERR "Keeping top $keep_top_P % of categories.\n";
}
elsif($keep_min_M){
	print STDERR "Removing low abundance and low count (<=$keep_min_M).\n";
}else{
	die "No keeping criteria specified.\n";
}

###############################################################################

open(FH, "<$input_st") || die "Could not open $input_st\n";

#------------------------------------------------------------------------------
# Read in header info
print STDERR "Loading header/category row...\n";
my $hdr_line=<FH>;
chomp $hdr_line;

my @hdr_col=split "\t", $hdr_line;

my @categories=@hdr_col;
shift @categories;	# Shift out sample_id
shift @categories;	# Shift out total

my $num_categories=$#categories+1;
print STDERR "Num categories found: $num_categories\n";
print STDERR "\n";

#------------------------------------------------------------------------------

my @sum_normalized;		# Normalized sums per category across all samples
my @sum_counts;			# Sum of counts per category across all samples
my %sample_counts_hash;		# Total counts per sample
my %all_sample_counts_hash;	# All count data

# Initialize counts, normalized, and totals
for(my $i=0; $i<$num_categories; $i++){
	$sum_normalized[$i]=0;
}

for(my $i=0; $i<$num_categories; $i++){
	$sum_counts[$i]=0;
}

#------------------------------------------------------------------------------

# Read in counts, but process 1 line at a time.
print STDERR "Loading samples...\n";
while(<FH>){
	chomp;

	my @data_col=split "\t", $_;

	# Remove first 2 columns
	my $sample_id=shift @data_col;
	my $total=shift @data_col;

	print STDERR "  $sample_id\n";

	my $sum=0;
	# Sum up counts
	for(my $i=0; $i<$num_categories; $i++){
		$sum+=$data_col[$i];
	}
	
	# Store sum per sample
	$sample_counts_hash{$sample_id}=$sum;

	# Normalize each category
	for(my $i=0; $i<$num_categories; $i++){

		my $normalized=$data_col[$i]/$sum;
		$sum_normalized[$i]+=$normalized;
		$sum_counts[$i]+=$data_col[$i];
	
	}
	
	# Store count data per sample
	$all_sample_counts_hash{$sample_id}=\@data_col;

}
print STDERR "done...\n";

###############################################################################

# Compute average across samples
my @samples=sort keys %sample_counts_hash;
my $num_samples=$#samples+1;

print STDERR "Num Samples Loaded: $num_samples\n";

# Compute average abundance across all samples
print STDERR "Computing average abundance per category...\n";
my %avgab_hash;
my %total_counts_hash;
for(my $i=0; $i<$num_categories; $i++){
	$avgab_hash{$categories[$i]}=$sum_normalized[$i]/$num_samples;
	$total_counts_hash{$categories[$i]}=$sum_counts[$i];
}

# Sort
print STDERR "Sorting by average abundance...\n";
my @sorted_categories = sort {$avgab_hash{$b} <=> $avgab_hash{$a}} @categories;

###############################################################################

# Determine keep criteria

my $keep_include_ix;

if($keep_top_N){
	$keep_include_ix=$keep_top_N-1;
}elsif($keep_top_P){
	print STDERR "Looking for top $keep_top_P% of categoies.\n";
	my $ab_cutoff=$keep_top_P/100.0;
	my $cumulative=0;
	for(my $i=0; $i<=$num_categories; $i++){
		my $cat=$sorted_categories[$i];
		$cumulative+=$avgab_hash{$cat};
		if($cumulative>=$ab_cutoff){
			$keep_include_ix=$i;
			print STDERR "Actual kept = " . ($cumulative*100) . " %\n";
			last;
		}
        }
}elsif($keep_min_M){
	# Search from the low end, for the first category that exceeds the min count.
	print STDERR "Looking for first non-singleton from end...\n";
	for(my $i=($num_categories-1); $i>=0; $i--){
		my $cat=$sorted_categories[$i];
		if($total_counts_hash{$cat}>$keep_min_M){
			$keep_include_ix=$i;
			last;
		}
	}
}

print STDERR "Top ", ($keep_include_ix + 1), " categories to be kept...\n";

###############################################################################

# Writing category/abundance/count file
if(0){
	print STDERR "Writing category/abundance/count file.\n";
	open(FH, ">$output_root\.abund.tsv") || die "Could not open $output_root\.abund.tsv\n";
	foreach my $cat(@sorted_categories){
		print FH "$cat\t$avgab_hash{$cat}\t$total_counts_hash{$cat}\n";
	}
	close(FH);
}

###############################################################################

# Look up has to get from sorted categories to counts.
my %cat_to_idx_hash;
my $i=0;
foreach my $cat(@categories){
	$cat_to_idx_hash{$cat}=$i;
	$i++;
}

# Output Summary Table 
open(FH, ">$output_root\.reduced.summary_table.tsv") || die "Could not open $output_root\.reduced.summary_table.tsv\n";

# Output Header Row
print FH "sample_id\ttotal";
for(my $i=0; $i<=$keep_include_ix; $i++){
	my $cat=$sorted_categories[$i];
	print FH "\t$cat";
}
print FH "\tRemaining\n";

# Output Samples
foreach my $samp_id(@samples){
	my $sample_tot=$sample_counts_hash{$samp_id};
	print FH "$samp_id\t$sample_tot";
	my $kept=0;
	for(my $i=0; $i<=$keep_include_ix; $i++){
		my $cat=$sorted_categories[$i];
		my $ix=$cat_to_idx_hash{$cat};
		my $count=${$all_sample_counts_hash{$samp_id}}[$ix];
		$kept+=$count;
		print FH "\t$count";
	}
	my $remaining=$sample_counts_hash{$samp_id}-$kept;
	print FH "\t$remaining\n";
}

###############################################################################

print STDERR "Done...\n";

