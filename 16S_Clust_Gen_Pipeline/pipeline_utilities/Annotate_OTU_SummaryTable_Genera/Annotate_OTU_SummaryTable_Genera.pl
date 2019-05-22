#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_s $opt_m $opt_o);

getopts("s:m:o:");
my $usage = "usage: 
$0 
	-s <Input OTU Summary Table>
	-m <OTU to Taxonomy Map, *.cons.taxonomy>
	-o <Output Genus Annotated OTU summary table>

	This script will read in the OTU summary table and the OTU to Taxa map,
	and output a new OTU table with the genera labeled.

";

if(!(
	defined($opt_s) && 
	defined($opt_m) && 
	defined($opt_o))){
	die $usage;
}

my $OTUTaxaMap=$opt_m;
my $SummaryTable=$opt_s;
my $STOut=$opt_o;

###############################################################################

print STDERR "Loading OTU to Taxonomy Map: $OTUTaxaMap\n";

my %otu_genus_map;

open(OTU_TAXA_MAP_FH, "<$OTUTaxaMap") || die "Could not open $OTUTaxaMap\n";

while(<OTU_TAXA_MAP_FH>){
	chomp;
	my ($otu_id, $size, $taxa)=split /\t/, $_;

	my @taxa_list=split ";", $taxa;
	my $genus=pop @taxa_list;

	$genus=~s/\(\d+\)$//;
	$genus=~s/\[/_/g;
	$genus=~s/\]/_/g;
	$genus=~s/\-//g;
	$genus=~s/_unclassified$/_uncl/;
	$genus=~s/_group$/_grp/;

	$otu_genus_map{$otu_id}="$otu_id.$genus";

	print STDERR "$otu_id / $genus: $otu_genus_map{$otu_id}\n";
	
}

close(OTU_TAXA_MAP_FH);

###############################################################################

print STDERR "Remapping OTU Summary Table: $SummaryTable\n";

open(ST_FH, "<$SummaryTable") || die "Could not open $SummaryTable\n";

open(ST_OUT_FH, ">$STOut") || die "Could not open $STOut\n";

my $header=1;
while(<ST_FH>){
	chomp;

	if($header){

		my $col_str=$_;
		my @cols=split "\t", $col_str;

		my @new_cols;

		# Copy over Sample_ID and Total
		push @new_cols, (shift @cols);
		push @new_cols, (shift @cols);

		# Remap IDs
		for(my $i=0; $i<=$#cols; $i++){
			my $orig=$cols[$i];
			my $trans=$otu_genus_map{$orig};
			if(!defined($trans)){
				print STDERR "Error, could not map OTU $orig to genus.\n";
				exit;
			}else{
				push @new_cols, $trans;
			}
		}
		
		# Output new header
		my $out_hdr=join "\t", @new_cols;
		print ST_OUT_FH "$out_hdr\n";
	
		$header=0;	
	}else{

		# Output counts, as is.
		print STDERR ".";
		print ST_OUT_FH "$_\n";	
	}

}

print STDERR "\n";

close(ST_FH);
close(ST_OUT_FH);

#------------------------------------------------------------------------------

print STDERR "Done.\n";
