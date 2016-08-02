#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_a $opt_b $opt_A $opt_B $opt_o);

getopts("A:B:a:b:o:");

my $usage = "
	usage:
	$0
		-a <Input filename, model A>
		-b <Input filename, model B>	
		-A <Stoichiometry Multiplier for A>
		-B <Stoichiometry Multiplier for B>
		-o <Output filename root>

	This script will read in two reaction files, identify
	the objective function, set their bounds to 0, and then
	create a new objective function that merges both
	objective functions with new stoichiometry based on
	specified multipliers.

";

if(!defined($opt_o) || 
	!defined($opt_a) || !defined($opt_A) ||
	!defined($opt_b) || !defined($opt_B)
){
	die $usage;
}

my $orgAfname=$opt_a;
my $orgBfname=$opt_b;
my $A_mult=$opt_A;
my $B_mult=$opt_B;
my $Output=$opt_o;

###############################################################################

print STDERR "  Organism A: $orgAfname\n";
print STDERR "Multiplier A: $A_mult\n";
print STDERR "  Organism B: $orgBfname\n";
print STDERR "Multiplier B: $B_mult\n";
print STDERR "Output: $Output\n";
print STDERR "\n";

my $FACTOR=1;
my $denom=$A_mult+$B_mult;
my $factorA=$FACTOR*$A_mult/$denom;
my $factorB=$FACTOR*$B_mult/$denom;

print STDERR "Factor A: $factorA\n";
print STDERR "Factor B: $factorB\n";

###############################################################################

my $OutputFn="$Output." . $A_mult . "v" . $B_mult . ".tsv";

open(OUT_FH, ">$OutputFn") || print STDERR "Could not open output file: $OutputFn.\n";

###############################################################################

my $RXN_ID_COL=0;
my $RXN_COL=2;
my $COMP_COL=4;
my $LB_COL=5;
my $UB_COL=6;
my $OBJ_COL=7;

my $RULE_COL=8;
my $SUBSYSTEM_COL=9;

my %rxn_id_hash;

sub ProcessReactions{
	my $filename=shift;
	my $keep_header=shift;

	open(IN_FH, "<$filename") || die "Could not open input file $filename\n";

	my $i=0;
	my $obj_rxn_str="";
	while(<IN_FH>){
		chomp;
		my @fields=split /\t/, $_;
		$#fields=9;

		if($fields[$OBJ_COL]==0){
			my $is_reaction=($fields[$RXN_COL]=~/>/);
			my $row_str=$_;
			if($is_reaction){
				if(defined($rxn_id_hash{$fields[$RXN_ID_COL]})){
					$fields[$RXN_ID_COL].="." . $fields[$COMP_COL];
					$row_str=join "\t", @fields;
				}
				$rxn_id_hash{$fields[$RXN_ID_COL]}=1;
				print OUT_FH "$row_str\n";
			}elsif($keep_header){
				print OUT_FH "$_\n";
			}
		}else{

			# Set current objective functions to off and bounds to 0
			$fields[$OBJ_COL]=0;
			$fields[$LB_COL]=0;
			$fields[$UB_COL]=0;
			#$fields[$RULE_COL]="";
			#$fields[$SUBSYSTEM_COL]="";
			my $outstr=join "\t", @fields;
			print OUT_FH "$outstr\n";

			$obj_rxn_str=$fields[$RXN_COL];
		}

		$i++;
	}
	print STDERR "Num Reactions Loaded: $i\n";
	close(IN_FH);

	return($obj_rxn_str);
}

###############################################################################

sub parse_additions{
	my $addition_str=shift;
	#print STDERR "$addition_str\n";
	my @metabolites=split /\+/, $addition_str;
	my @parsed_coef;
	my @parsed_meta;
	foreach my $metabolite(@metabolites){
		$metabolite=~s/^ //g;
		$metabolite=~s/ $//g;
		#print "$metabolite\n";

		my $coef=1;
		my $meta="";
		if($metabolite=~/\(([0-9\.]+)\)\s+(.+)/){
			$coef=$1;
			$meta=$2;
		}else{
			$meta=$metabolite;	
		}
	
		#print STDERR "$coef / $meta\n";
		push @parsed_coef, $coef;
		push @parsed_meta, $meta;	
	}
	return(\@parsed_coef, \@parsed_meta);
}

sub parse_reaction{
	my $reaction=shift;	

	my ($compart, $reaction)=split " : ", $reaction;
	
	my $compart_p;
	if($compart=~/\[(\S+)\]/){
		$compart_p=$1;
		print STDERR "Compartment: $compart_p\n";
	}else{
		print STDERR "Error: Could not parse out compartment.\n";
		print STDERR "Reaction: '$reaction'\n";
		die;
	}

	my ($lhs, $rhs)=split "-->", $reaction;
	#print "$lhs\n";
	my ($lhs_coef, $lhs_meta)=parse_additions($lhs);
	#print "$rhs\n";
	my ($rhs_coef, $rhs_meta)=parse_additions($rhs);

	return($lhs_coef, $lhs_meta, $rhs_coef, $rhs_meta, $compart_p);
}

sub multi_coef{
	my $arr_ref=shift;
	my $val=shift;

	for(my $i=0; $i<=$#{$arr_ref}; $i++){
		${$arr_ref}[$i]*=$val;
	}
}

sub to_string{
	my $coef_ref=shift;
	my $meta_ref=shift;
	my $compart=shift;

	my @coef=@{$coef_ref};
	my @met=@{$meta_ref};
	
	my @combo;

	for(my $i=0; $i<=$#coef; $i++){
		#print STDERR "-----------------\n";
		#print STDERR "$coef[$i]\n";
		#print STDERR "$met[$i]\n";
		#print STDERR "$compart\n";
		#print STDERR "-----------------\n";
		if($coef[$i]!=0){
			push @combo, ("($coef[$i]) $met[$i]" . "[$compart]");	
		}
	}

	my $str=join " + ", @combo;
	return($str);
}

sub combine{
	my $A_lhs_coef=shift; my $A_lhs_meta=shift; my $A_rhs_coef=shift; my $A_rhs_meta=shift; my $A_compart=shift; my $A_mult=shift;
	my $B_lhs_coef=shift; my $B_lhs_meta=shift; my $B_rhs_coef=shift; my $B_rhs_meta=shift; my $B_compart=shift; my $B_mult=shift;

	multi_coef($A_lhs_coef, $A_mult);
	multi_coef($A_rhs_coef, $A_mult); 
	multi_coef($B_lhs_coef, $B_mult);
	multi_coef($B_rhs_coef, $B_mult); 

	my $A_lhs=to_string($A_lhs_coef, $A_lhs_meta, $A_compart);
	my $A_rhs=to_string($A_rhs_coef, $A_rhs_meta, $A_compart);
	my $B_lhs=to_string($B_lhs_coef, $B_lhs_meta, $B_compart);
	my $B_rhs=to_string($B_rhs_coef, $B_rhs_meta, $B_compart);

	my $lhs;
	my $rhs;

	if($A_lhs ne ""  && $B_lhs ne ""){
		$lhs="$A_lhs + $B_lhs";
	}
	if($A_lhs eq ""  && $B_lhs ne ""){
		$lhs="$B_lhs";
	}
	if($A_lhs ne ""  && $B_lhs eq ""){
		$lhs="$A_lhs";
	}

	if($A_rhs ne ""  && $B_rhs ne ""){
		$rhs="$A_rhs + $B_rhs";
	}
	if($A_rhs eq ""  && $B_rhs ne ""){
		$rhs="$B_rhs";
	}
	if($A_rhs ne ""  && $B_rhs eq ""){
		$rhs="$A_rhs";
	}

	my $complete_str="$lhs --> $rhs";

	return($complete_str);

}

###############################################################################

my $A_reaction=ProcessReactions($orgAfname, 1);
my $B_reaction=ProcessReactions($orgBfname, 0);

# combineRXN CombinedObjectiveFunction Irreversible CombinedCompart 0 1000 1 
my ($A_lhs_coef, $A_lhs_meta, $A_rhs_coef, $A_rhs_meta, $A_compart)= parse_reaction($A_reaction);
my ($B_lhs_coef, $B_lhs_meta, $B_rhs_coef, $B_rhs_meta, $B_compart)= parse_reaction($B_reaction);

my $combined_rxn=combine($A_lhs_coef, $A_lhs_meta, $A_rhs_coef, $A_rhs_meta, $A_compart, $factorA,
			$B_lhs_coef, $B_lhs_meta, $B_rhs_coef, $B_rhs_meta, $B_compart, $factorB);

print STDERR "\nReaction Equation:\n";
print STDERR "$combined_rxn\n";

my $reaction_str=join "\t", (
	"Combined_RXN",				# Reaction ID
	"Combined_Objective_$A_mult:$B_mult",	# Reaction Description
	$combined_rxn,				# Stoichiometry
	"Irreversible",				# Reversibility
	"$A_compart,$B_compart",		# Compartment
	"0",					# LowerBound
	"1000", 				# UpperBound
	"1",					# Objective?
	"",					# Gene rules
	""					# Subsystem
);

print OUT_FH "$reaction_str\n";

close(OUT_FH);

###############################################################################

print STDERR "\n";
print STDERR "Done.\n";
