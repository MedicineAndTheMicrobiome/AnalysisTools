#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
#use lib "$FindBin::Bin";
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_s $opt_u $opt_a $opt_o $opt_v $opt_b);

getopts("s:u:a:o:v:b:");

my $UB_CUTOFF=.8;
my $AB_CUTOFF=.01;

my $MINOR_UB_CUTOFF=.5;
my $MINOR_AB_CUTOFF=.01;

my $library_path="$FindBin::Bin/lib";

my $usage = "usage:
$0
	-s <Sample summary_table.tsv>
	[-u <Ubiquity cutoff, default=$UB_CUTOFF>]
	[-a <Abundance cutoff, default=$AB_CUTOFF>]

	[-v <Minor Core Ubiquity cutoff, default=$MINOR_UB_CUTOFF>]
	[-b <Minor Core Abundance cutoff, default=$MINOR_AB_CUTOFF>]

	[-o <output filename root>]

	This script runs all the individual sample Corbata analyses.
	These include:
		a.) Compute number of core at specified Ubiquity/Abundance
		b.) Make core Ub-Ab plots
			color based on core
			remaining is grey
		c.) Make minor core Ub-Ab plot
		d.) Plot abundance vs. variance plots

	Dependendent libraries should be at:
		$library_path
        
";

if(!(
	defined($opt_s) 
)){
	die $usage;
}

###############################################################################

my $ubiquity_cutoff=$UB_CUTOFF;
my $abundance_cutoff=$AB_CUTOFF;

my $minor_core_ubiquity_cutoff=$MINOR_UB_CUTOFF;
my $minor_core_abundance_cutoff=$MINOR_AB_CUTOFF;

my $output_root;

if(defined($opt_u)){
	$ubiquity_cutoff=$opt_u;
}
if(defined($opt_a)){
	$abundance_cutoff=$opt_a;
}
if(defined($opt_o)){
	$output_root=$opt_o;
}
if(defined($opt_v)){
	$minor_core_ubiquity_cutoff=$opt_v;
}
if(defined($opt_b)){
	$minor_core_abundance_cutoff=$opt_b;
}


my $sample_summary_table=$opt_s;

if(!defined($output_root)){
	my ($name, $path)=fileparse($sample_summary_table);
	
	$name=~s/\.summary_table\.tsv$//;
	$name=~s/\.summary_table\.xls$//;
	$name=~s/\.summary_table\.txt$//;

	$output_root=$name;
}

print STDERR "Input Sample: $sample_summary_table\n";
print STDERR "Output Root: $output_root\n";

###############################################################################

sub execute{
	my $cmd=shift;
	$cmd=~s/\s+/ /g;
	print STDERR "Executing $cmd\n";
	my $result=`$cmd`;
	print STDERR $result;
}

###############################################################################

sub compute_core{
	my $summary_table=shift;
	my $ubiquity=shift;
	my $abundance=shift;
	my $output_root=shift;

	my $cmd="
		$library_path/Compute_Number_of_Core.r
			-i $summary_table
			-u $ubiquity
			-a $abundance
			-o $output_root
	";
	execute("$cmd");
	return;
};


sub generate_color_map{
	my $category_list=shift;

	my $out_name=$category_list;
	$out_name=~s/\.txt$//;
	$out_name=~s/\.tsv$//;

	my $cmd="
		$library_path/Assign_Colors.r
			-i $category_list
			-o $out_name
	";
	execute("$cmd");
	return;
};

sub generate_variation_plot{
	my $st=shift;
	my $output_root=shift;

	my $cmd="
		$library_path/TaxonomicVariation.r
			-i $st
			-o $output_root
	";
	execute("$cmd");
	return;
};

sub generate_UbAb_plot{
	my $st=shift;
	my $ub_cutoff=shift;
	my $ab_cutoff=shift;
	my $minor_ub_cutoff=shift;
	my $minor_ab_cutoff=shift;
	my $output=shift;
	my $smooth=shift;	# Smoothing: 0 off, 1 on
	my $idmap=shift;
	my $color_sch=shift;
	my $keep_list=shift;

	# Generate CDFs
	my $options="";
	if($smooth){
		$options.=" -s";
	}

	print STDERR "$st\n";
	my $cmd="
		$library_path/Compute_CDFs.r 
			-i $st
			-o $output
			$options
	";
	execute("$cmd");	

	# Generate name to differentiate smooth/observed
	my $cdf_file=$output;
	$cdf_file.=".lin.log";
	if($smooth){
		$cdf_file.=".smoothed"
	}
	$cdf_file.=".cdf";

	# Plot major core
	my $cmd="
		$library_path/Plot_CDFs.r
			-i $cdf_file
			-u $ub_cutoff
			-a $ab_cutoff
			-c $color_sch
			-r
	";
	execute("$cmd");
	
	# Plot minor and other core
	my $cmd="
		$library_path/Plot_CDFs_OtherCore.r
			-i $cdf_file
			-u $minor_ub_cutoff
			-a $minor_ab_cutoff
			-c $color_sch
			-r
	";
	execute("$cmd");

}


###############################################################################

compute_core($sample_summary_table, $ubiquity_cutoff, $abundance_cutoff, $output_root);
# Generates:
# 	$output_root.core_statistics.tsv
# 	$output_root.core_members.tsv
#	$output_root.observed_core.tsv
#	$output_root.significant_core.tsv


generate_color_map(
	"$output_root.significant_core.tsv"
);
# 	$output_root.significant_core.colormap
# 	$output_root.significant_core.colormap.pdf

generate_color_map(
	"$output_root.observed_core.tsv"
);
# 	$output_root.observed_core.colormap
# 	$output_root.observed_core.colormap.pdf


# Smoothed UbAb Plot
generate_UbAb_plot(
	$sample_summary_table,
	$ubiquity_cutoff,
	$abundance_cutoff,
	$minor_core_ubiquity_cutoff,
	$minor_core_abundance_cutoff,
	"$output_root",
	1,
	"",
	"$output_root.observed_core.colormap",
	"$output_root.observed_core.tsv"
);

# Observed UbAb Plot
generate_UbAb_plot(
	$sample_summary_table,
	$ubiquity_cutoff,
	$minor_core_abundance_cutoff,
	$minor_core_ubiquity_cutoff,
	$abundance_cutoff,
	"$output_root",
	0,
	"",
	"$output_root.observed_core.colormap",
	"$output_root.observed_core.tsv"
);

# Generate variation vs abundance plot
generate_variation_plot(
	$sample_summary_table,
	$output_root
);

print STDERR "Done.\n";
