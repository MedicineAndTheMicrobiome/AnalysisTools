#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_n $opt_o $opt_c);

getopts("n:o:c:");
my $usage = "usage: 

$0 

	-n <count_table file>
	-o <output file>
	<-c \"Comment\">


	This script will read in a .count_table file from mothur
	and give the number of representatives and the numbers
	represented.
	
	The results will be appended to the specified output file name.
	
	The -c option can be used to describe the name of the result
	appended. If the -c option is not specified then the input
	file will be used.

";

if(!(
	defined($opt_o) &&
	defined($opt_n))){
	die $usage;
}

my $count_table_file=$opt_n;
my $output_file=$opt_o;
my $comment=$opt_c;

if(!defined($comment)){
	my ($fname, $path)=fileparse($names_file);
	$comment=$fname;
}

print STDERR "Count Table File: $count_table_file\n";
print STDERR "Output File: $output_file\n";

###############################################################################

my $print_hdr=1;
if(-e $output_file){
	$print_hdr=0;	
}

open(FH, "<$count_table_file") || die "Could not open $count_table_file.\n";

 my $num_representatives=0;
 my $num_represented=0;

# To remove first 3 lines of the file
my $skiphead = 3;
<FH> for (1..$skiphead);

# To count number of representative and number of represented from the count_table
while (<FH>) {

    my ($representatives, $represented, $countmatrix) = split;

    $num_representatives += 1;
    $num_represented += $represented;

}

# print "The number of represented is $num_represented\n";
# print "The number of representative is $num_representatives\n";

close(FH);

my $collapse_rate=$num_representatives/$num_represented;

###############################################################################

open(OUT_FH, ">>$output_file") || die "Could not open $output_file.\n";

if($print_hdr){
	print OUT_FH "Comment\tNum_Reptvs\tNum_Reptd\tCol_Rate\n";
}

print OUT_FH "$comment\t$num_representatives\t$num_represented\t$collapse_rate\n";

close(OUT_FH);

###############################################################################

print STDERR "done.\n";
