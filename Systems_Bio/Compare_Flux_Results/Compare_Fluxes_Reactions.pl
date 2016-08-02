#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_a $opt_b $opt_o);

getopts("a:b:o:");

my $usage = "
	usage:
	$0
		-a <Input flux filename, Primary>
		-b <Input flux filename, Secondary>
		-o <Output comparison filename root>

	This script will compare the results of two FBA results.

	The columns should be:
		ID			Reaction ID/Abbreviation
		Description		Reaction Description/Name
		Flux value		Output from FBA
		Compartment
		Equation		Implied direction (If Flux Value is <0, then equation should go from Left <-- Right

	The Fluxes will be with respect to the Primary flux.  This means
	that all primary fluxes will be positive.  Reactions will always
	go from left to right.  
	
";

if(!defined($opt_a) || !defined($opt_b) || !defined($opt_b)){
	die $usage;
}

my $Primary=$opt_a;
my $Secondary=$opt_b;
my $Output=$opt_o;

###############################################################################

my %compartment_map;
$compartment_map{"c0"}="c";
$compartment_map{"e0"}="e";

###############################################################################

print STDERR "Primary: $Primary\n";
print STDERR "Secondary: $Secondary\n";
print STDERR "Output: $Output\n";

###############################################################################

sub reaction_to_string{
	my $reaction_str=shift;
	my $compartment=shift;

	# Clean spaces
	print "Reaction: '$reaction_str'\n";
	$reaction_str=~s/\(\)//g;
	$reaction_str=~s/  / /g;
	$reaction_str=~s/^ +//;
	$reaction_str=~s/ +$//;
	print "Clean Reaction: '$reaction_str'\n";
	
	# Identify reversibility
	my $reversible;
	my $split_str;
	my $swapped=0;

	if($reaction_str=~/-->/){
		$split_str="-->";
		$swapped=0;
	}elsif($reaction_str=~/<--/){
		$split_str="<--";
		$swapped=1;
	}else{
		die "Could not detect reversibility in reaction: $reaction_str\n";
	}

	# Identify left and right hand side
	$reaction_str=~s/\s//g;
	my ($lhs, $rhs)=split $split_str, $reaction_str;

	# Swapp lhs with rhs, so irreversible reaction is always going to the right
	if($swapped){
		($lhs, $rhs)=($rhs, $lhs)
	}

	# Append global compartment to end of each metabolite
	my @lhs_arr=sort(split "\\+", $lhs);
	my @rhs_arr=sort(split "\\+", $rhs);

	if($compartment ne ""){
		# Translate global compartment i.e., [c0] : ...  into [c] : ..
		if(defined($compartment_map{$compartment})){
			$compartment=$compartment_map{$compartment};
			print STDERR "$compartment\n";
		}
		for(my $i=0; $i<=$#lhs_arr; $i++){
			$lhs_arr[$i].=$compartment;
		}
		for(my $i=0; $i<=$#rhs_arr; $i++){
			$rhs_arr[$i].=$compartment;
		}
	}else{
		# Translate individual compartments i.e., [c0] into [c];
		my $ind_comp;
		for(my $i=0; $i<=$#lhs_arr; $i++){
			$lhs_arr[$i]=~/\[(.+)\]/;
			$ind_comp=$compartment_map{$1};
			if(defined($ind_comp)){
				$lhs_arr[$i]=~s/\[(.+)\]/[$ind_comp]/;
			}
		}
		for(my $i=0; $i<=$#rhs_arr; $i++){
			$rhs_arr[$i]=~/\[(.+)\]/;
			$ind_comp=$compartment_map{$1};
			if(defined($ind_comp)){
				$rhs_arr[$i]=~s/\[(.+)\]/[$ind_comp]/;
			}
		}
	}

	# Reform string
	$lhs=join "+", @lhs_arr;
	$rhs=join "+", @rhs_arr;

	print "$lhs-->$rhs\n\n";

	return($lhs, $rhs);

}

###############################################################################

my $ABB_COL=0;
my $NAM_COL=1;
my $FLUX_COL=2;
my $COMP_COL=3;
my $EQN_COL=4;

my %reaction_hash;

sub load_reaction{
	my $filename=shift;

	print STDERR "Loading: $filename\n";
	open(IN_FH, "<$filename") || die "Could not open input file $filename\n";

	my %reaction_hash;
	my %flux_hash;
	
	my $i=0;
	while(<IN_FH>){
		chomp;

		my $reaction;
		my $compartment="";
		
		my @fields=split "\t", $_;

		my $rxn_abbreviation=$fields[$ABB_COL];
		my $rxn_name=$fields[$NAM_COL];
		my $rxn_flux=$fields[$FLUX_COL];
		my $rxn_comp=$fields[$COMP_COL];
		my $rxn_eqn=$fields[$EQN_COL];

		# Skip if row is the header 
		if(lc($rxn_eqn) eq "equation"){
			next;
		}
		#print "$rxn_field\n";

		# Parse out the compartment info
		my $compartment;
		if($rxn_eqn=~/ : /){
			($compartment, $reaction)=split " : ", $rxn_eqn;
		}else{
			($compartment, $reaction)=("", $rxn_eqn);
		}
		
		# Standardize the LHS, RHS 
		my ($lhs, $rhs)=reaction_to_string($reaction, $compartment);

		# Store both direction so we can differentiate between mismatching direction and non-existance
		$reaction_hash{"$lhs-->$rhs"}="$rxn_abbreviation\t1";
		$reaction_hash{"$rhs-->$lhs"}="$rxn_abbreviation\t-1";

		# Keep the flux
		$rxn_flux=abs($rxn_flux);	# The negative pos doesn't matter since we are using the arrows for direction
		$flux_hash{$rxn_abbreviation}="$rxn_flux\t$lhs-->$rhs";

		# Keep track of the number of reactions loaded
		$i++;
	}

	print STDERR "Num Reactions Loaded: $i\n";
	close(IN_FH);

	return(\%reaction_hash, \%flux_hash);
}

###############################################################################

my ($prim_reaction_hash, $prim_flux_hash)=load_reaction($Primary);
my ($sec_reaction_hash, $sec_flux_hash)=load_reaction($Secondary);

#--------------------------------------------------------------------------------

my @prim_rxn_id_arr=sort keys %{$prim_flux_hash};
my @sec_rxn_id_arr=sort keys %{$sec_flux_hash};

my @matching;
my @primary_only;
my @secondary_only;

my %secondary_hits;
my @primary_misses;

foreach my $pri_rxn_id (@prim_rxn_id_arr){

	my $primary_flux_info=${$prim_flux_hash}{$pri_rxn_id};
	my ($pri_flux, $pri_eqn)=split "\t", $primary_flux_info;
	my $sec_rxn_info=${$sec_reaction_hash}{$pri_eqn};

	if(defined($sec_rxn_info)){

		my ($sec_rxn_id, $direction)=split "\t", $sec_rxn_info;
		my $secondary_flux_info=${$sec_flux_hash}{$sec_rxn_id};
		my ($sec_flux, $sec_eqn)=split "\t", $secondary_flux_info;

		$sec_flux*=$direction;

		my $match_info="$pri_rxn_id\t$pri_flux\t$sec_rxn_id\t$sec_flux\t$pri_eqn";
		print STDERR "$match_info\n";
		push @matching, $match_info;

		# Keep track of which secondary id were hit
		$secondary_hits{$sec_rxn_id}=1;

	}else{
		# Keep track of which primary id didn't hit anything
		push @primary_misses, $pri_rxn_id;
	}
}

###############################################################################

close(REPORT_FH);
close(ACTION_FH);

###############################################################################

my $primary_only="$Output\.primary_only.tsv";
my $secondary_only="$Output\.secondary_only.tsv";
my $shared="$Output\.shared.tsv";

#-----------------------------------------------------------------------------

open(PRI_ONLY_FH, ">$primary_only") || die "Could not open $primary_only\n";
foreach my $id(@primary_misses){
	print PRI_ONLY_FH "$id\t${$prim_flux_hash}{$id}\n";
}
close(PRI_ONLY_FH);

#-----------------------------------------------------------------------------

open(SEC_ONLY_FH, ">$secondary_only") || die "Could not open $secondary_only\n";
foreach my $id (@sec_rxn_id_arr){
	if(!defined($secondary_hits{$id})){
		print SEC_ONLY_FH "$id\t${$sec_flux_hash}{$id}\n";
	}
}
close(SEC_ONLY_FH);

#-----------------------------------------------------------------------------

open(SHARED_FH, ">$shared") || die "Could not open $shared\n";
foreach my $rxn_info(sort @matching){
	print SHARED_FH "$rxn_info\n";

}

#-----------------------------------------------------------------------------

close(PRI_ONLY_FH);
close(SEC_ONLY_FH);
close(SHARED_FH);

###############################################################################


print STDERR "\n";
print STDERR "Done.\n";
