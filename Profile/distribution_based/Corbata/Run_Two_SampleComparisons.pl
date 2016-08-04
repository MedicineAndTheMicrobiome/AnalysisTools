#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
#use lib "$FindBin::Bin";
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_a $opt_b $opt_o);

getopts("a:b:o:");

my $library_path="$FindBin::Bin/lib";

my $usage = "usage:
$0
	-a <Sample A summary_table.tsv>
	-b <Sample B summary_table.tsv>
	[-o <output filename root>]

	This script will perform all the two sample comparisons available
	in the Corbata suite of tools:
		1.) Generate taxa specific comparison statistics
			a.) Identify categories with significance <.05 and <.1
		2.) Generate U-U Plot
			a.) Compute CDFs (for each sample)
			b.) Compare/Draw UU Plots (smoothed and unsmoothed)
				i.)  Only drawing taxa with significance <.1
				ii.) Only coloring taxa with significance <.05
		3.) Compute AWKS Statistic
			a.) generates p-values and null distribution plots

	Dependendent libraries should be at:
		$library_path
        
";

if(!(
	defined($opt_a) && 
	defined($opt_b) 
)){
	die $usage;
}

###############################################################################

my $sample_a_summary_table=$opt_a;
my $sample_b_summary_table=$opt_b;
my $output_root=$opt_o;

if(!defined($output_root)){
	my ($aname, $apath)=fileparse($sample_a_summary_table);
	my ($bname, $bpath)=fileparse($sample_b_summary_table);
	
	$aname=~s/\.summary_table\.tsv$//;
	$bname=~s/\.summary_table\.tsv$//;

	$aname=~s/\.summary_table\.xls$//;
	$bname=~s/\.summary_table\.xls$//;

	$aname=~s/\.summary_table\.txt$//;
	$bname=~s/\.summary_table\.txt$//;

	$output_root="$aname\_vs_$bname";
}

print STDERR "Input Sample A: $sample_a_summary_table\n";
print STDERR "Input Sample B: $sample_b_summary_table\n";
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

sub compare_taxa{
	my $sample_a_summary_table=shift;
	my $sample_b_summary_table=shift;
	my $output_root=shift;

	my $cmd="
		$library_path/Compare_Microbiomes.r
			-a $sample_a_summary_table
			-b $sample_b_summary_table
			-o $output_root
	";
	execute("$cmd");
	return;
};

sub compute_AWKS{
        my $sample_a_summary_table=shift;
        my $sample_b_summary_table=shift;
        my $output_root=shift;

	my $cmd="
		$library_path/Compare_wAWKS.r
			-a $sample_a_summary_table
			-b $sample_b_summary_table
			-o $output_root
	";
	execute("$cmd");
	return;
}

sub generate_color_map{
	my $category_list=shift;

	my $out_name=$category_list;
	$out_name=~s/\.txt$//;

	my $cmd="
		$library_path/Assign_Colors.r
			-i $category_list
			-o $out_name
	";
	execute("$cmd");
	return;
};

sub generate_UU_plot{
	my $st_a=shift;
	my $st_b=shift;
	my $output=shift;
	my $smooth=shift;	# Smoothing: 0 off, 1 on
	my $idmap=shift;
	my $color_sch=shift;
	my $keep_list=shift;
	my $fn_ext=shift;

	#----------------------------------------------------------------------
	# Generate CDFs
	my @st_files=($st_a, $st_b);
	my @cdf_files;

	my $options="";
	if($smooth){
		$options.=" -s";
	}

	foreach my $st(@st_files){
		print STDERR "$st\n";
		my $cmd="
			$library_path/Compute_CDFs.r 
				-i $st
				$options
		";
		execute("$cmd");	

		my $cdf_file=$st;
		$cdf_file=~s/\.summary_table\.tsv$//;
		$cdf_file=~s/\.summary_table\.xls$//;

		$cdf_file.=".lin.log";
		if($smooth){
			$cdf_file.=".smoothed"
		}
		$cdf_file.=".cdf";

		push @cdf_files, $cdf_file;
	}
	

	#----------------------------------------------------------------------
	# Generate CDFs
	# Generate Plots
	my $options="";
	if($idmap ne ""){
		$options.=" -m $idmap";
	}
	if(defined($idmap)){
		$options.=" -c $color_sch";
	}
	if($keep_list ne ""){
		$options.=" -k $keep_list";
	}

	$output.=".$fn_ext";

	my $cmd="
		$library_path/Compare_CDFs.r
			-A $cdf_files[0]
			-B $cdf_files[1]
			-o $output
			$options
	";
	
	execute("$cmd");

	return;
}

###############################################################################

compare_taxa(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root
);
# Generates:
# 	$output_root.compare.pdf
# 	$output_root.exc_dif.txt
# 	$output_root.pdf
# 	$output_root.p-values.txt
# 	$output_root.ubiq_diff.txt
#	$output_root.uncor_pval_lt05
#	$output_root.uncor_pval_lt10
# 	$output_root.weighted_ks.txt

generate_color_map(
	"$output_root.uncor_pval_lt05"
);
# Generate:
# 	$output_root.uncor_pval_lt05.colormap.pdf
# 	$output_root.uncor_pval_lt05.colormap

#----------------------------------------------------------------

generate_UU_plot(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root,
	0,
	"",
	"$output_root.uncor_pval_lt05.colormap",
	"$output_root.uncor_pval_lt10",
	""
);
# Generates:
# 	$output_root.UU.pdf

generate_UU_plot(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root,
	1,
	"",
	"$output_root.uncor_pval_lt05.colormap",
	"$output_root.uncor_pval_lt10",
	"smoothed"
);
# Generates:
# 	$output_root.smoothed.UU.pdf

#----------------------------------------------------------------

generate_UU_plot(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root,
	0,
	"",
	"$output_root.uncor_pval_lt05.colormap",
	"",
	"grey"
);
# Generates:
# 	$output_root.UU.grey.pdf

generate_UU_plot(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root,
	1,
	"",
	"$output_root.uncor_pval_lt05.colormap",
	"",
	"smoothed.grey"
);
# Generates:
# 	$output_root.smoothed.UU.grey.pdf

#----------------------------------------------------------------
	
compute_AWKS(
	$sample_a_summary_table,
	$sample_b_summary_table,
	$output_root
);
# Generates:
# 	$output_root.awks_test.pdf
# 	$output_root.awks_test.tsv

print STDERR "Done.\n";
