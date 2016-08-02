#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o $opt_m $opt_c $opt_n $opt_u);

getopts("i:o:c:m:n:u:");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-o <Translated output filename>
		-m <Map file, tab-separated>
		-c <compartment name to translate, for example: 'e0'>
		-n <new compartment name>
		[-u <used ids>]

	This script will read in a reactions file, and then translate all
	the compounds that are in the specified compartment.

	The map file should have the format:
		<key>\\t<new value>\\n

";

if(!defined($opt_i) || !defined($opt_o) || !defined($opt_m) || !defined($opt_c) || !defined($opt_n)){
	die $usage;
}

my $InputFile=$opt_i;
my $OutputFile=$opt_o;
my $MapFile=$opt_m;
my $Compartment=$opt_c;
my $NewCompartmentName=$opt_n;
my $OutputUsedIDs=$opt_u;

###############################################################################

print STDERR "Input: $InputFile\n";
print STDERR "Map: $MapFile\n";
print STDERR "Compartment: $Compartment\n";

###############################################################################

my %map;
my $num_map_entries=0;
open(MAP_FH, "<$MapFile") || die "Could not open $MapFile\n";
while(<MAP_FH>){
	chomp;
	my ($key, $val)=split /\t/, $_;

	if($val ne ""){
		$map{$key}=$val;
	}else{
		$map{$key}="UNMAPPED";
	}
	$num_map_entries++;
}
close(MAP_FH);

print STDERR "Num Map Entries Read: $num_map_entries\n";

###############################################################################

delete($map{""});
my @keys=keys %map;

my %used_ids;

sub translate{
	my $in=shift;

	print STDERR "BEFORE: $in\n";

	# Single compartment reaction
	if($in=~/^\[(\S+)\]/){
		if($1 eq $Compartment){
			#print STDERR "Whole Reaction\n";
			$in=~s/^\[\S+\]/\[$NewCompartmentName\]/;
			my @parts=split / /, $in;
			for(my $i=0; $i<=$#parts; $i++){
				my $new_id=$map{$parts[$i]};
				if(defined($new_id)){
					$parts[$i]=$new_id;
					$used_ids{$new_id}=$NewCompartmentName;
				}
			}
			$in=join " ", @parts;
		}else{
			#print STDERR "Other Compartment.\n";
		}
	# Multiple compartment reaction
	}else{
		#print STDERR "Individual Compounds\n";
		while($in=~/(\S+)\[$Compartment\]/){
			my $old_id=$1;
			my $new_id=$map{$old_id};
			$in=~s/$old_id\[$Compartment\]/$new_id\[$NewCompartmentName\]/;
			$used_ids{$new_id}=$NewCompartmentName;
			#print STDERR "$in\n";
		}

	}
	print STDERR "AFTER: $in\n\n";

	return($in);
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";
open(OUT_FH, ">$OutputFile") || die "Could not open input file $OutputFile\n";


my $fh;
my $temp_filename;

my $REACTION_COLUMN=2;
my $COMPARTMENT_COLUMN=4;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;

	# Grab reaction column
	my $reaction=$array[$REACTION_COLUMN];
	# Remap reaction compounds and location
	my $new_reaction=translate($reaction);
	# Reinsert reaction back into array
	$array[$REACTION_COLUMN]=$new_reaction;


	# Rename compartment
	my $compartments_str=$array[$COMPARTMENT_COLUMN];
	$compartments_str=~s/\s+//g;
	$compartments_str=~s/"//g;
	my @comparts=split /,/, $compartments_str;
	for(my $i=0; $i<=$#comparts; $i++){
		if($comparts[$i] eq $Compartment){
			$comparts[$i]=$NewCompartmentName;
		}
	}
	$compartments_str=join ",", @comparts;
	$array[$COMPARTMENT_COLUMN]=$compartments_str;

	# Output new reaction 
	$#array=9;
	my $outstr=join "\t", @array;	
	print OUT_FH "$outstr\n";

	# Output pulse
	$input_lines++;
	if(!($input_lines%10)){
		print STDERR ".";
	}
}

close(IN_FH);
close(OUT_FH);

###############################################################################

if(defined($OutputUsedIDs)){
	open(USED_FH, ">$OutputUsedIDs") || die "Could not open $OutputUsedIDs for writing.\n";
	foreach my $new_id(keys %used_ids){
		print USED_FH "$new_id\t$used_ids{$new_id}\n";
	}
	close(USED_FH);
}

###############################################################################

print STDERR "\n";
print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Done.\n";
