#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_t);

getopts("t:");

my $usage = "
	usage:
	$0
		-t <TREMBL dat file>

	Reads in UniProt TREMBL dat file and generates a simplified
	table of ID, description, PFAM ID, TIGRFAM ID, Taxa ID, GO IDs, and length.

	Output goes to STDOUT.
	
";

if(!defined($opt_t)){
	die $usage;
}

my $tremble_dat_filename=$opt_t;

###############################################################################

my $cur_AC="";
my $cur_OC="";
my $cur_ID="";
my @pfam;
my @tigrfam;
my @desc;
my $taxa_id="";
my @ec;
my @go_p;
my @go_f;
my $length;

open(FH, "<$tremble_dat_filename") || die "Could not open $tremble_dat_filename\n";


my $i=0;
while(<FH>){
	chomp $_;
	
	if($_=~/^AC   /){
		# Extract accession/UniRef ID
		my ($tag, $rec)=split "   ", $_;
		$rec=~s/;$//;
		$cur_AC=$rec;

	}elsif($_=~/^DE   /){

		# Extract description and EC if available
		if($_=~/SubName: Full=(.+);/){
			push @desc, $1;
		}elsif($_=~/RecName: Full=(.+);/){
			push @desc, $1;
		}elsif($_=~/EC=(.+);/){
			push @ec, $1;
		}

	}elsif($_=~/^DR   /){
	
		# Extract pfam, tiger fram, and go ids
		my ($tag, $rec)=split "   ", $_;
		my @subrec=split "; ", $rec;

		if($subrec[0] eq "Pfam"){
			push @pfam, $subrec[1];	
		}elsif($subrec[0] eq "TIGRFAMs"){
			push @tigrfam, $subrec[1];	
		}elsif($subrec[0] eq "GO"){
			my $gotype=substr($subrec[2],0,1);
			if($gotype eq "P"){	
				push @go_p, $subrec[1];
			}elsif($gotype eq "F"){
				push @go_f, $subrec[1];
			}else{
				if($gotype ne "C"){
					print STDERR "GO parsing error: '$_'\n";
				}
			}
		}

	}elsif($_=~/^OX   /){

		# Extract taxa ids
		my ($tag, $rec)=split "   ", $_;
		my @subrec=split "; ", $rec;
		
		if($subrec[0]=~/NCBI_TaxID=(\d+)/){
			$taxa_id=$1;
		}

	}elsif($_=~/^SQ   /){
		
		# Extract sequence lengths
		if($_=~/SEQUENCE\s+(\d+) AA;/){
			$length=$1;
		}

	}elsif($_=~/^\/\//){

		# Dumping record out.

		# Convert arrays into strings
		my $pfamlist=join ";", sort(@pfam);
		my $tigrfamlist=join ";", sort(@tigrfam);
		my $delist=join ";", @desc;
		my $goplist=join ";", sort(@go_p);
		my $goflist=join ";", sort(@go_f);
		my $eclist=join ";", sort(@ec);

		my $outstr=join "\t", ($cur_AC, $delist, $taxa_id, $length, $pfamlist, $tigrfamlist, $goplist, $goflist, $eclist);

		print "$outstr\n";

		$length=$taxa_id=$cur_AC=$cur_OC="";
		@pfam=();
		@tigrfam=();
		@desc=();
		@go_p=();
		@go_f=();
		@ec=();

		# Heart beat
		if(!($i%250000)){
			print STDERR "$i Records Processed.\n";
		}
		$i++;
	}
}
close(FH);

###############################################################################

print STDERR "Done.\n";
