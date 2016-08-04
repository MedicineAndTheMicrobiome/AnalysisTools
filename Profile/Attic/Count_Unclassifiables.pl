#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_p $opt_d);
getopts("p:d:");


my $usage = "usage:
$0
	-p <project>
	[-d <alternate directory to write files>]

	reads in files:

		<project>.domain.summary_table.xls
		<project>.phylum.summary_table.xls
		<project>.class.summary_table.xls
		<project>.order.summary_table.xls
		<project>.family.summary_table.xls
		<project>.genus.summary_table.xls

	Generates a table for each sample of the number of reads that could not be classified.

	Writes out:

		<project>.level_counts.xls
		<project>.level_percentages.xls


";

if(!(
	defined($opt_p)
)){
	die $usage;
}

my $project_id=$opt_p;
my $outdir;

if(defined($opt_d)){
	$outdir=$opt_d;
}else{
	$outdir="";
}

my @levels=("domain", "phylum", "class", "order", "family", "genus");
my %count_hash;
my %sample_hash;

foreach my $level(@levels){

	my $input_file_name="$project_id\.$level\.summary_table.xls";

	open(INFH, "<$input_file_name") || die "Could not open $input_file_name\n";

	my $first_line=<INFH>;
	$first_line=~s/\n$//;
	my @categories=split(/\t/, $first_line);

	my $unknown_col=-1;
	my $col=0;

	print "Counting on $level\n";
	while(<INFH>){

		chomp;
		my @splitline=split /\t/,$_;

		my $sample_id=$splitline[0];
		my $sample_count=$splitline[1];
		
		if($level eq "domain"){
			$count_hash{"all"}{$sample_id}=$sample_count;
		}
		$count_hash{$level}{$sample_id}=$sample_count;
		$sample_hash{$sample_id}=1;
	}

	close(INFH);

}


#######################################################################

my $num_samples=keys %sample_hash;

my $count_fn="$outdir/$project_id\.level_counts.xls";
open(CNT_FN, ">$count_fn") || die "Could not open $count_fn for writing.\n";

print CNT_FN "Sample\t";
print CNT_FN join ("\t", ("all",@levels));
print CNT_FN "\n";

my %weighted_sum_hash;
foreach my $sample(sort keys %sample_hash){
	print CNT_FN "$sample";
	foreach my $level(("all",@levels)){
		print CNT_FN "\t$count_hash{$level}{$sample}";
		$weighted_sum_hash{$level}+=$count_hash{$level}{$sample};
	}
	print CNT_FN "\n";

}

print CNT_FN "[Nominal Average]";
foreach my $level(("all", @levels)){
	printf CNT_FN "\t%.2f", 1.0*$weighted_sum_hash{$level}/$num_samples;
}
print CNT_FN "\n";

print CNT_FN "[Weighted Average Percentage]";
foreach my $level(("all", @levels)){
	printf CNT_FN "\t%.2f", 100.0*$weighted_sum_hash{$level}/($weighted_sum_hash{"all"});
}
print CNT_FN "\n";




close(CNT_FN);

#######################################################################

my $perc_fn="$outdir/$project_id\.level_percentages.xls";
open(PRC_FN, ">$perc_fn") || die "Could not open $perc_fn for writing.\n";

print PRC_FN "Sample\t";
print PRC_FN join ("\t", ("all",@levels));
print PRC_FN "\n";

my %unweighted_sum_hash;
foreach my $sample(sort keys %sample_hash){
	print PRC_FN "$sample";
	foreach my $level(("all",@levels)){
		if($count_hash{"all"}{$sample} == 0){
			print STDERR "ignored: $sample\n";
		}else{
			my $perc=100*$count_hash{$level}{$sample}/$count_hash{"all"}{$sample};
			printf PRC_FN "\t%3.2f", $perc;
			$unweighted_sum_hash{$level}+=$perc;
		}
	}
	print PRC_FN "\n";

}

print PRC_FN "[Unweighted Average Percentage]";
foreach my $level(("all", @levels)){
	printf PRC_FN "\t%.2f", 1.0*$unweighted_sum_hash{$level}/$num_samples;
}
print PRC_FN "\n";

close(PRC_FN);

#######################################################################

print STDERR "Done.\n\n";
