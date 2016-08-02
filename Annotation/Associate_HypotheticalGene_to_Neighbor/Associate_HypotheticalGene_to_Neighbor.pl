#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_f);

getopts("f:o:");

my $usage = "
	usage:
	$0
		-f <GenBank Feature Table>

	Reads in GenBank Feature table file and extract description for
	each locus id and generates a map.  If a protein is hypothetical, then 
	it reports the description of the closest up and downstream gene.

	Output goes to STDOUT.
	
";

if(!defined($opt_f)){
	die $usage;
}

my $feature_table=$opt_f;

###############################################################################

my @records;

sub to_record{
	my $begin=shift;
	my $end=shift;
	my $locus_id=shift;
	my $gene=shift;
	my $gene_syn=shift;	
	my $product=shift;
	my $EC=shift;

	my $EC_str="";
	if($EC ne ""){
		$EC_str=" (EC:$EC)";
	}

	my $gene_str="";
	if($gene ne ""){
		$gene_str.="$gene;";
	}
	if($gene_syn ne ""){
		$gene_str.="$gene;";
	}
	if($gene_str ne ""){
		$gene_str.=" ";
	}
	
	my $info="$gene_str$product$EC_str";
	if($info eq ""){
		$info="[No Description]";
	}

	if($begin>$end){
		($begin, $end)=($end, $begin);
	}
	my $bps=($end-$begin)+1;

	my $info_string="$begin-$end\t$locus_id\t$info {$bps bp}";			
	push @records, $info_string;
}

###############################################################################

open(FH, "<$feature_table") || die "Could not open $feature_table\n";

my $begin="";
my $end="";

my $prev_begin="";
my $prev_end="";

my $locus_id;
my $gene;
my $product;
my $EC;
my $gene_syn;

while(<FH>){
	chomp;

	if($_=~/^>/){
		next;
	}else{
		my @cols=split "\t", $_;
		if(($cols[0] eq "") && ($cols[1] eq "") && ($cols[2] eq "")){
			my $tag=$cols[3];
			my $info=$cols[4];
	
			if($tag eq "locus_tag"){
				$locus_id=$info;
			}elsif($tag eq "gene"){
				$gene=$info;
			}elsif($tag eq "product"){
				$product=$info;
			}elsif($tag eq "EC_number"){
				$EC=$info;
			}
				
		}else{
			$begin=$cols[0];
			$end=$cols[1];

			#if(($cols[2] ne "gene") && 
			#	($cols[2] ne "CDS") && 
			#	($cols[2] ne "rRNA") && 
			#	($cols[2] ne "tRNA") && 
			#	($cols[2] ne "tmRNA") && 
			#	($cols[2] ne "misc_RNA") 
			#){
			#	print STDERR "Third column in position record not gene or CDS tag: '$cols[2]'\n";
			#}

			if(($prev_begin ne $begin || $prev_end ne $end) &&($prev_end ne "" && $prev_begin ne "")){
				to_record($begin,$end,$locus_id,$gene,$gene_syn,$product,$EC);			

				$locus_id=$gene=$product=$EC=$gene_syn="";
	
			}

			$prev_begin=$begin;
			$prev_end=$end;
		}
			

	}
}
to_record($begin,$end,$locus_id,$gene,$gene_syn,$product,$EC);			

close(FH);

###############################################################################
# Sort by the begin coordinate

sub by_coords{
	my ($a_begin, $a_rest)=split "-", $a;
	my ($b_begin, $b_rest)=split "-", $b;
	return($a_begin <=> $b_begin);
}

my @sorted_records=sort by_coords @records;

###############################################################################
# Go through each record and identity hypothetical genes and associate them with neighbors

my $num_records=$#sorted_records;
my @hypothetical_annotated;

for(my $i=0; $i<$num_records; $i++){
	my $rec=$sorted_records[$i];
	my ($coord, $loc_id, $info)=split "\t", $rec;

	#my $info_string="$begin-$end\t$locus_id\t$info {$bps bp}";			

	if($info=~/^hypothetical protein/){

		# Grab up and downstream gene info
		my ($coord_up, $loc_id_up, $info_up)=split "\t", $sorted_records[$i-1];
		my ($coord_down, $loc_id_down, $info_down)=split "\t", $sorted_records[$i+1];

		# Substitute description with locus id
		$info=~s/hypothetical protein/$loc_id/;
		$info_up=~s/hypothetical protein/$loc_id_up/;
		$info_down=~s/hypothetical protein/$loc_id_down/;

		# 
		$info="$info [-$info_up / +$info_down]";
	}

	$hypothetical_annotated[$i]="$loc_id\t$info";
}

###############################################################################
# Output record to stdout

foreach my $rec (@hypothetical_annotated){
	print STDOUT "$rec\n";
}

###############################################################################

print STDERR "Done.\n";
