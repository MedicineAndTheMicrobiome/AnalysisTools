#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o $opt_m $opt_c $opt_s);

getopts("i:o:c:m:s");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-o <Translated output filename>
		-m <Map file, tab-separated>
		[-s (swap mapping so second column is the key)]
		-c <list of Columns to translate, starting from 0>

	This script will parse the reaction formula and translate
	the ID with exact match according to the map file.

	The map file should have the format:
		<key>\\t<new value>\\n

";

if(!defined($opt_i) || !defined($opt_m) || !defined($opt_c)){
	die $usage;
}

my $InputFile=$opt_i;
my $OutputFile=$opt_o;
my $MapFile=$opt_m;
my $ColumnNum=$opt_c;
my $SwapMapKey=defined($opt_s);

###############################################################################

print STDERR "\n";
print STDERR "Input: $InputFile\n";
print STDERR "Map: $MapFile\n";
print STDERR "Column Number: $ColumnNum\n";
print STDERR "Swap Key/Map Columns: $SwapMapKey\n";
print STDERR "\n";

###############################################################################

my %map;
my %degenerate_hash;
my $num_map_entries=0;

open(MAP_FH, "<$MapFile") || die "Could not open $MapFile\n";
while(<MAP_FH>){
	chomp;
	my ($key, $val)=split /\t/, $_;

	if($SwapMapKey){
		($key, $val)=($val, $key);
	}

	if($val ne ""){
		if(!defined($map{$key})){
			$map{$key}=$val;
		}else{
			print STDERR "Duplicate key/val mapping: $key -> $val (Also mapped to: $map{$key})\n";
			$degenerate_hash{$key}=1;
		}
	}else{
		$map{$key}="UNMAPPED";
	}
	$num_map_entries++;
}
close(MAP_FH);

print STDERR "Num Map Entries Read: $num_map_entries\n";
print STDERR "\n";

###############################################################################

delete($map{""});

sub parse_equation{
	my $eq=shift;
	my $dir;
	
	if($eq=~/<==>/){
		$dir="<==>";	
	}elsif($eq=~/-->/){
		$dir="-->";
	}elsif($eq=~/<--/){
		$dir="<--";
	}

	my ($lhs, $rhs)=split / $dir */, $eq;
	return($lhs, $rhs, $dir);	
}	

sub parse_stoichiometric_components{
	my $stcomp=shift;
	my @ctcomps=split / \+ /, $stcomp;
	return(\@ctcomps);
}

sub construct_stoichiometric_components{
	my $comp_ref=shift;
	my $str=join " + ", @{$comp_ref};
	return($str);
}

sub parse_indiv_stoich_comp{
	my $ind_stcomp=shift;
	my ($coef, $met, $compart);
	if($ind_stcomp=~/\(([0-9\.]+)\) (.+)\[(.+)\]$/){
		($coef, $met, $compart)=($1, $2, $3);
	}elsif($ind_stcomp=~/(.+)\[(.+)\]$/){
		($coef, $met, $compart)=(1, $1, $2);
	}elsif($ind_stcomp=~/\(([0-9\.]+)\) (.+)/){
		($coef, $met, $compart)=($1, $2, "");
	}elsif($ind_stcomp=~/(.+)/){
		($coef, $met, $compart)=(1, $1, "");
	}else{
		print STDERR "Error parsing $ind_stcomp.\n";
	}

	#print STDERR "Before: $ind_stcomp Interpreted As: ($coef) $met [$compart]\n";

	return($coef, $met, $compart);

}

sub construct_indiv_stoich_comp{
	my $coef=shift;
	my $met=shift;
	my $compart=shift;

	my $str="";
	if($coef!=1){
		$str.="($coef) ";
	}
	$str.=$met;
	if($compart ne ""){
		$str.="[$compart]";	
	}
	return($str);
}

sub translate_metabolite{
	my $met=shift;

	if(defined($degenerate_hash{$met})){
		die "Error! Requested to translate degenerate metabolite: $met\n";
	}

	if(defined($map{$met})){
		$met=$map{$met};
	}else{
		print STDERR "Untranslated: $met\n";
	}
	return($met);
}

sub translate_reaction{
	my $reaction_field=shift;

	my ($left, $right)=split / : /, $reaction_field;
	my $global_comp;
	my $equation;
	if($right eq ""){
		$global_comp="";
		$equation=$left;
	}else{
		$global_comp=$left;
		$equation=$right;
	}

	my ($lhs, $rhs, $dir)=parse_equation($equation);
	my $lhs_stoich_comp=parse_stoichiometric_components($lhs);
	my $rhs_stoich_comp=parse_stoichiometric_components($rhs);

	for(my $i=0; $i<=$#{$lhs_stoich_comp}; $i++){
		my ($coef, $met, $compart)=parse_indiv_stoich_comp(${$lhs_stoich_comp}[$i]);
		$met=translate_metabolite($met);
		${$lhs_stoich_comp}[$i]=construct_indiv_stoich_comp($coef, $met, $compart);
	}

	for(my $i=0; $i<=$#{$rhs_stoich_comp}; $i++){
		my ($coef, $met, $compart)=parse_indiv_stoich_comp(${$rhs_stoich_comp}[$i]);
		$met=translate_metabolite($met);
		${$rhs_stoich_comp}[$i]=construct_indiv_stoich_comp($coef, $met, $compart);
	}

	my $new_lhs=construct_stoichiometric_components($lhs_stoich_comp);
	my $new_rhs=construct_stoichiometric_components($rhs_stoich_comp);

	my $new_reaction="$new_lhs $dir $new_rhs";

	if($global_comp ne ""){
		$new_reaction="$global_comp : $new_reaction";
	}	
	return($new_reaction);
	
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";
open(OUT_FH, ">$OutputFile") || die "Could not open input file $OutputFile\n";


my $fh;
my $temp_filename;
my $max_col=0;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;

	my $num_col=$#array+1;
	$max_col=($max_col>$num_col)?$max_col:$num_col;

	my $str_to_translate=$array[$ColumnNum];
	my $translated_str=translate_reaction($str_to_translate);

	#print STDERR "\n";
	#print STDERR "$str_to_translate\n";
	#print STDERR "$translated_str\n";

	$array[$ColumnNum]=$translated_str;
	my $joined=join "\t", @array;
	$joined.= "\t" x ($max_col-$num_col);
	print OUT_FH "$joined\n";

	$input_lines++;
}

close(IN_FH);
close(OUT_FH);

###############################################################################

print STDERR "\n";
print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Done.\n";
