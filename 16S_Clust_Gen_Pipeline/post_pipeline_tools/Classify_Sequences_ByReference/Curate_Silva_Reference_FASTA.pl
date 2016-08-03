#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_g $opt_o);

getopts("f:g:o:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-g <target genus, e.g. Enterococcus>
	-o <Output Filename Root>

	This script will extract all the sequences that look like
	they have been properly assigned to a species for the targeted genus.

	The output is a reduced fasta file and a mapping file.

	Records with the following assignments will be removed:
		<genus>_sp.* 
		_bacterium*
		uncultured*
		unidentified*
	
	Strain info will be removed.
		e.g. Enterococcus faecalis EnGen0215 will be changed to Enterococcus faecalis

	The map will contain an accession to genus_species map.

	The resultant FASTA file can be used as a reference database for BLAST.
	This fasta will only contain the kept sequences and only the accession will
	be in the defline.
	

";

if(!(
	defined($opt_f) && 
	defined($opt_g) && 
	defined($opt_o))){
	die $usage;
}


my $input_fasta=$opt_f;
my $target_genus=$opt_g;
my $output_fn_root=$opt_o;

###############################################################################
# Read in features

my $num_found=0;

open(FASTA_FH, "<$input_fasta") || die "Could not open $input_fasta\n";
open(FASTA_OUT_FH, ">$output_fn_root\.kept.fasta") || die "Could not open $output_fn_root\.kept.fasta";
open(KEEP_OUT_FH, ">$output_fn_root\.kept.map") || die "Could not open $output_fn_root\.kept.map";
open(REMOVE_OUT_FH, ">$output_fn_root\.rem.map") || die "Could not open $output_fn_root\.rem.map";

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

close(FASTA_FH);
close(FASTA_OUT_FH);
close(KEEP_OUT_FH);
close(REMOVE_OUT_FH);

print STDERR "$num_found kept\n";
print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $subgenus;
	my $genus_species;
	my $skip=0;
	my $accession;
	if($defline=~/^>(\S+)\s+(.*)/){
		$accession=$1;
		my $description=$2;

		my @taxa=split ";", $description;
		$subgenus=$taxa[6];

		if($subgenus=~/^uncultured/ || $subgenus=~/metagenome/){
			$skip=1;	
		}elsif($subgenus=~/^unidentified/){
			$skip=1;
		}elsif($subgenus=~/^bacterium/ || $subgenus=~/\s+bacterium/){
			$skip=1;
		}elsif($subgenus=~/$target_genus sp\./){
			$skip=1;
		}else{
			# remove strain info:
			if($subgenus=~/$target_genus (\S+)/){
				$genus_species="$target_genus $1";
			}else{
				$skip=1;
			}
		}
	}else{
		die "Error parsing defline: $defline\n";
	}

	if($skip==1){
		print STDOUT "Removing: $defline\n";
		print REMOVE_OUT_FH "$accession\t$defline\n";
	}else{
		print STDOUT "Keeping: $subgenus -> $genus_species\n";
		print KEEP_OUT_FH "$accession\t$genus_species\n";

		# Output FASTA record
		print FASTA_OUT_FH ">$accession\n";
		my $length=length($sequence);
		my $width=50;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print FASTA_OUT_FH substr($sequence, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);
		
		$num_found++;
	}
}

#------------------------------------------------------------------------------

