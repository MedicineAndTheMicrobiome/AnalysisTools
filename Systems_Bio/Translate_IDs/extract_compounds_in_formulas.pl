#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c);

getopts("i:c:");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-c <Columns starting from 0 which contains the formulas>

	This script will pull out all the compounds in the reaction equations
	and report them in STDOUT.

";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $InputFile=$opt_i;
my $ColumnNum=$opt_c;

###############################################################################

print STDERR "\n";
print STDERR "Input: $InputFile\n";
print STDERR "Column Number: $ColumnNum\n";
print STDERR "\n";

###############################################################################

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

my %compounds_found;

sub pull_compounds_from_reaction{
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
		$compounds_found{$met}=1;
	}

	for(my $i=0; $i<=$#{$rhs_stoich_comp}; $i++){
		my ($coef, $met, $compart)=parse_indiv_stoich_comp(${$rhs_stoich_comp}[$i]);
		$compounds_found{$met}=1;
	}
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

while(<IN_FH>){
	chomp;

	my @array=split /\t/, $_;
	my $rxn_formula=$array[$ColumnNum];
	pull_compounds_from_reaction($rxn_formula);

	$input_lines++;
}

close(IN_FH);

###############################################################################

foreach my $cmpd (sort keys %compounds_found){
	print STDOUT "$cmpd\n";
}

###############################################################################

print STDERR "\n";
print STDERR "Done.\n";
