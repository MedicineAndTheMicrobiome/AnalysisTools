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
		-a <Input filename, Primary>
		-b <Input filename, Secondary>
		-o <Output filename root>

	This script will read in two reaction files which
	have already been mapped to use the same IDs (i.e. BiGG),
	and then based on the input and output stoichiometries
	tell you if they are the same and if they have the same
	reversibilities.

";

if(!defined($opt_a) || !defined($opt_b) || !defined($opt_o)){
	die $usage;
}

my $Primary=$opt_a;
my $Secondary=$opt_b;
my $Output=$opt_o;

###############################################################################

print STDERR "Primary: $Primary\n";
print STDERR "Secondary: $Secondary\n";
print STDERR "Output: $Output\n";

###############################################################################

sub reaction_to_string{
	my $reaction_str=shift;
	my $compartment=shift;

	#print "Reaction: '$reaction_str'\n";
	
	# Identify reversibility
	my $reversible;
	my $split_str;
	if($reaction_str=~/-->/){
		$reversible=0;
		$split_str="-->";
	}elsif($reaction_str=~/<==>/){
		$reversible=1;
		$split_str="<==>";
	}else{
		die "Could not detect reversibility in reaction: $reaction_str\n";
	}

	# Identify left and right hand side
	$reaction_str=~s/\s//g;
	my ($lhs, $rhs)=split $split_str, $reaction_str;

	# Append compartment
	my @lhs_arr=sort(split "\\+", $lhs);
	my @rhs_arr=sort(split "\\+", $rhs);

	for(my $i=0; $i<=$#lhs_arr; $i++){
		$lhs_arr[$i].=$compartment;
	}
	for(my $i=0; $i<=$#rhs_arr; $i++){
		$rhs_arr[$i].=$compartment;
	}
	
	# Reform string
	$lhs=join "+", @lhs_arr;
	$rhs=join "+", @rhs_arr;

	#print "$lhs\n";
	#print "$rhs\n";
	#print "$reversible\n";	

	return($lhs, $rhs, $reversible);

}

###############################################################################

my $ABB_COL=0;
my $NAM_COL=1;
my $RXN_COL=2;
my $REV_COL=3;

my %reaction_hash;

sub load_reaction{
	my $filename=shift;

	open(IN_FH, "<$filename") || die "Could not open input file $filename\n";

	my @clean_reactions;
	my $i=0;
	while(<IN_FH>){
		chomp;

		my $reaction;
		my $compartment="";
		
		my @fields=split "\t", $_;

		my $rxn_abbreviation=$fields[$ABB_COL];
		my $rxn_name=$fields[$NAM_COL];
		my $rxn_field=$fields[$RXN_COL];

		# Skip if row is the header 
		if($rxn_field eq "equation"){
			next;
		}
		print "$rxn_field\n";

		# Parse out the compartment info
		my $compartment;
		if($rxn_field=~/ : /){
			($compartment, $reaction)=split " : ", $rxn_field;
		}else{
			($compartment, $reaction)=("", $rxn_field);
		}
		
		
		# Standardize the LHS, RHS and decode reversibility
		my ($lhs, $rhs, $reversible)=reaction_to_string($reaction, $compartment);

		# Keep the reaction information together
		push @clean_reactions, (join " / ", ($rxn_abbreviation, $rxn_name, $lhs, $rhs, $reversible));

		# Keep track of the number of reactions loaded
		$i++;
	}

	print STDERR "Num Reactions Loaded: $i\n";
	close(IN_FH);

	return(\@clean_reactions);
}

###############################################################################

my $prim_arr_ref=load_reaction($Primary);
my $sec_arr_ref=load_reaction($Secondary);

# Compare

#--------------------------------------------------------------------------------

my $out_fname="$Output\.rxn_match\.tsv";
open(REPORT_FH, ">$out_fname") || die "Could not open $out_fname\n";

my $hdr_str=join "\t", (
	"ReferenceModel",
	"SecondaryModel",
	"SecondaryRxnAbbr",
	"MatchingDirectionality",
	"RefReversibility",
	"SecReversibility",
	"MatchingReversibilities",
	"NeedsCuration",
	"MakeSecRev",
	"MakeSecIrr",
	"RevSecMakeIrr",
	"Contradiction"
);
print REPORT_FH "$hdr_str\n";

#--------------------------------------------------------------------------------

my $out_action_fname="$Output\.action\.tsv";
open(ACTION_FH, ">$out_action_fname") || die "Could not open $out_action_fname\n";

$hdr_str=join "\t", (
	"SecondaryRxnAbbr",
	"ReferenceRxnAbbr",
	"TargetDirection",
	"TargetReversibility"
);
print ACTION_FH "$hdr_str\n";

#--------------------------------------------------------------------------------

# For each reaction in the reference model...
foreach my $pri_str (@{$prim_arr_ref}){

	my ($rxn_abbr1, $nam1, $lhs1, $rhs1, $rev1)=split " / ", $pri_str;

	# Look through each reaction in the secondary model for a match
	foreach my $sec_str (@{$sec_arr_ref}){

		my ($rxn_abbr2, $nam2, $lhs2, $rhs2, $rev2)=split " / ", $sec_str;

		my $match=0;
		my $direction_match=0;
		
		# Determine if there are any reactions that are the same only based on stoichiometry
		if(($lhs1 eq $lhs2) && ($rhs1 eq $rhs2)){
			#print "$nam1 == $nam2	(Same Direction)\t(Same Reversibility)\n";
			$match=1;
			$direction_match=1;
		}elsif(($lhs1 eq $rhs2) && ($rhs1 eq $lhs2)){
			#print "$nam1 == $nam2	(Reverse Direction)\n";
			$match=1;
			$direction_match=-1;
		}
		
		# Output a line for each reaction match
		if($match){
		
			# Do reversibilities match?
			my $matching_rev=0;
			if($rev1==$rev2){
				$matching_rev=1;	
			}

			#--------------------------------------------------------------------------
			# Flag the reactions that are ok, so we can identify reactions that are not
			# ok, but can't explain
		
			# Assume reaction needs curation unless...
			my $needs_curation=1;
			#    Direction and reversibilities match
			if($direction_match==1 && $matching_rev){
				$needs_curation=0;
			}	
			# or Directions don't match but both are reversible
			if($direction_match==-1 && ($rev1==1 && $rev2==1)){
				$needs_curation=0;
			}

			#--------------------------------------------------------------------------
			# These fixes require changing reversibility only 

			# Determine action to get secondary model to match reference model
			my $make_sec_rev=0;
			my $make_sec_irr=0;

			# If both reactions going in the same direction, and m1 is not reversible, but m2 is reversible
			if($direction_match==1 && $rev1==0 && $rev2==1){
				# Strengthen secondary model
				$make_sec_irr=1;
			}

			# If m1 is reversible but m2 is not reversible
			if($rev1==1 && $rev2==0){
				# Loosen secondary model
				$make_sec_rev=1;
			}

			#--------------------------------------------------------------------------
			# These fixes require changing reaction direction and reversibility

			# If m1 and m2 are going opposite directions, but m1 demands irreversibility
			my $rev_sec_make_irr=0;
			if($direction_match==-1 && ($rev1==0) && ($rev2==1)){
				# Swap reaction directionality and make irreversible
				$rev_sec_make_irr=1;
			}

			# If m1 and m2 are going opposite direction, and each demands irreversibility
			my $contradiction=0;
			if($direction_match==-1 && ($rev1==0) && ($rev2==0)){
				# Flag as contraction
				$contradiction=1;
			}

			#--------------------------------------------------------------------------
			# Output a report

			my $out_str=join "\t", (
				$nam1,
				$nam2,
				$rxn_abbr2,
				$direction_match,
				$rev1,
				$rev2,
				$matching_rev,
				$needs_curation,
				$make_sec_rev,
				$make_sec_irr,
				$rev_sec_make_irr,
				$contradiction
			);

			print REPORT_FH "$out_str\n";

			if($needs_curation){
				my $out_str=join "\t", (
					$rxn_abbr2,
					$rxn_abbr1,
					$direction_match,
					$rev1				
				);

				print ACTION_FH "$out_str\n";
			}
	
		}
			

	}
}



#foreach my $str (@{$prim_arr_ref}){
#	print "$str\n";
#}

#foreach my $str (@{$sec_arr_ref}){
#	print "$str\n";
#}

###############################################################################

close(REPORT_FH);
close(ACTION_FH);

print STDERR "\n";
print STDERR "Done.\n";
