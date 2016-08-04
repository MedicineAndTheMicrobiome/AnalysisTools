#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;

use vars qw($opt_i $opt_o);
getopts("i:o:");

my $usage="
	
	$0 
		-i <input bergeyTrainingTree.xml, eg \$SVN/16sDataAnalysis/trunk/site_analysis/aux/java/data/classifier/bergeyTrainingTree.xml>
		-o <bergey taxon names/ranks>

	This script will read in the bergeyTrainingTree and generate a list of all of 
	the ranks that could be generated from the from the xml file.

	Output will look like:
		<rank>\\t<name>\\n
		<rank>\\t<name>\\n
		<rank>\\t<name>\\n
		...
		<rank>\\t<name>\\n

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;	
}

my $xml=$opt_i;
my $outroot=$opt_o;

###############################################################################

open(FH, "<$xml") || die "Could not open $xml for reading.\n";

my %tree;
my %rank_hash;

while(<FH>){
	
	chomp;
	$_=~s/^\<TreeNode //;
	$_=~s/\>\<\/TreeNode\>$//;

	if($_=~/name="(.+)" taxid="(\d+)" rank="(\S+)" parentTaxid="(\d+)" leaveCount="(\d+)" genusIndex="([\d-]+)"/){
		my ($name, $taxid, $rank, $parent_taxid, $leave_count, $genus_index)=($1, $2, $3, $4, $5, $6);
		#print "$name, $taxid, $rank, $parent_taxid, $leave_count, $genus_index)\n";
		$name=~s/&quot;/"/g;
		$name=~s/ /_/g;
		@{$tree{$taxid}}=($name, $rank, $parent_taxid);
		push @{$rank_hash{$rank}}, $taxid;
	}else{
		if(!(($_=~/trainsetNo/) || ($_=~/name="Root"/))){
			print "$_\n";
			print "\tparse error\n";
		}
	}
}	

###############################################################################

print STDERR "Ranks found:\n";
foreach my $rank(sort keys %rank_hash){
	print STDERR "\t$rank\n"
}
print STDERR "\n";

###############################################################################

my @ranks_of_interest=(keys %rank_hash);

open(OUT_ALL_FH, ">$outroot\.all") || die "Could not open $outroot\.all\n";
open(OUT_NOSUB_FH, ">$outroot\.no_sub") || die "Could not open $outroot\.no_sub\n";

foreach my $rank(@ranks_of_interest){

	print STDERR "Working on $rank\n";

	foreach my $member(@{$rank_hash{$rank}}){
		
		my ($full_name, $subfree_name)=climb_to_root($member);
		print OUT_ALL_FH "$rank\t$full_name\n";

		if($full_name ne $subfree_name){
			print OUT_ALL_FH "$rank\t$subfree_name\n";
		}

		if(!($rank=~/^sub/)){
			print OUT_NOSUB_FH "$rank\t$subfree_name\n";
		}
	}

}

close(OUT_ALL_FH);
close(OUT_NOSUB_FH);

print STDERR "\ndone.\n\n";

###############################################################################

sub climb_to_root{
	my $taxid=shift;
	
	my @full_name;
	my @subfree_name;

	my ($name, $rank, $parent_taxid);
	do{
		($name, $rank, $parent_taxid)=@{$tree{$taxid}};

		unshift @full_name, $name;
		if(!($rank=~/^sub/)){
			unshift @subfree_name, $name;
		}

		$taxid=$parent_taxid;
	}while($rank ne "domain");
	
	my $full_name=join " ", @full_name;
	my $subfree_name=join " ", @subfree_name;
	return($full_name, $subfree_name);
}

