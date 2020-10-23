#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

my $FILES_TO_COLUMN_BIN="~/git/AnalysisTools/Column/files_to_columns.pl";

getopts("i:");

my $usage = "
	usage:
	$0

	-i <output directory of where all block runs>

";

if(
	!defined($opt_i)
){
	die $usage;
}


my $OutputDir=$opt_i;
my $OutputSummaryDir="$OutputDir\.summaries";

###############################################################################

my @extensions=(
	"div_as_pred.anova.summary.tsv",
	"div_as_resp.anova.summary.tsv",
	"alr_as_pred.anova.summary.tsv",
	"alr_as_resp.anova.summary.tsv",
	"perm.anova.summary.tsv");


###############################################################################

print STDERR "\n";
print STDERR "Output Directory:          $OutputDir\n";
print STDERR "Output Summary Directory:  $OutputSummaryDir\n";
print STDERR "\n";

###############################################################################

if(! -e $OutputSummaryDir){
	mkdir $OutputSummaryDir;
}

my @files=split /\n/, `find $OutputDir | grep summary.tsv`;

foreach my $ext (@extensions){
	my @targets;

	print "Look for $ext files...\n";
	foreach my $file(@files){
		if($file=~/$ext$/){
			push @targets, $file;
		}
	}
	print "\n";

	print "Accumulating...\n";
	for my $targ(@targets){
		print "\t$targ\n";
	}
	my $targ_list=join " ", @targets;


	`$FILES_TO_COLUMN_BIN $targ_list > $OutputSummaryDir/$ext`;

	print "\n";

}

