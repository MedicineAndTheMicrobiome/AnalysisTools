#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_a $opt_o);

getopts("i:a:o:");

my $usage = "
	usage:
	$0
		-i <Input Model>
		-a <Action File>
		-o <Output Model>

	This script will apply the modification in the action file
	to the input model.  

";

if(!defined($opt_i) || !defined($opt_a) || !defined($opt_o)){
	die $usage;
}

my $Model=$opt_i;
my $Action=$opt_a;
my $Output=$opt_o;

###############################################################################

print STDERR "Model: $Model\n";
print STDERR "Action: $Action\n";
print STDERR "Output: $Output\n";

###############################################################################

sub load_actions{
	my $fname=shift;

	my %action_hash;

	open(FH, "<$fname") || die "Could not open $fname.\n";

	while(<FH>){
		chomp;
		my ($mod_abbr, $ref_abbr, $tar_dir, $tar_rev)=split "\t", $_;	
		my $rec=join "\t", ($ref_abbr, $tar_dir, $tar_rev);
		$action_hash{$mod_abbr}=$rec;
	}		

	close(FH);

	return(\%action_hash);
}

###############################################################################


my $ABB_COL=0;
my $NAM_COL=1;
my $EQU_COL=2;
my $REV_COL=3;
my $COM_COL=4;
my $LB_COL =5;
my $UB_COL =6;
my $OBJ_COL=7;
my $RUL_COL=8;
my $SUB_COL=9;

###############################################################################

my $action_hash_ref=load_actions($Action);

open(IN_FH, "<$Model") || die "Could not open $Model for reading.\n";
open(OUT_FH, ">$Output") || die "Could not open $Output for writing.\n";

while(<IN_FH>){
	chomp;
	my @fields=split "\t", $_;
	my $out_str;

	if($fields[$EQU_COL] eq "equation"){
		print OUT_FH "$_\n";
	}else{
		my $react_abb=$fields[$ABB_COL];
		my $react_rec=${$action_hash_ref}{$react_abb};
		
		# If this reaction needs to be updated ...
		if(defined($react_rec)){

			my ($ref_abbr, $tar_dir, $tar_rev)=split "\t", $react_rec;	
			my @new_fields=@fields;

			# Identify reaction dir/reversibility
			my $rxn_eqn=$fields[$EQU_COL];
			my $dir;
			if($rxn_eqn=~/-->/){
				$dir=" --> ";	
			}elsif($rxn_eqn=~/<==>/){
				$dir=" <==> ";
			}
			
			# Parse out only the reaction part of the string by removing
			# the compartment information
			my $single_comp=0;
			my $clean_eqn=$rxn_eqn;
			my $compartment;
			if($rxn_eqn=~/ : /){
				$single_comp=1;	
				($compartment, $clean_eqn)=split / : /, $rxn_eqn;
			}

			# Determine left and right hand side of the reaction
			my ($lhs, $rhs)=split /$dir/, $clean_eqn;	

			# Determine symbol for target reversiblity/irreversibility
			my $new_dir;
			if($tar_rev == 0){
				$new_dir=" --> ";
			}elsif($tar_rev == 1){
				$new_dir=" <==> ";
			}

			# If we need to swap the sides, do so now
			if($tar_dir == -1){
				($lhs, $rhs)=($rhs, $lhs);
			}

			# Genereate new reaction string
			my $new_rxn_eqn=$lhs . $new_dir . $rhs;

			# Reconstruct equation with compartment
			if($single_comp){
				$new_rxn_eqn=$compartment . " : " . $new_rxn_eqn;
			}
			
			# Set reversibility
			my $rev_str;
			my $lb_str;
			my $ub_str;
			if($tar_rev==0){
				$rev_str="Irreversible";
				$lb_str=0;
				$ub_str=1000;
			}else{
				$rev_str="Reversible";
				$lb_str=-1000;
				$ub_str=1000;
			}

			# Package up new fields
			$new_fields[$LB_COL]=$lb_str;
			$new_fields[$UB_COL]=$ub_str;
			$new_fields[$REV_COL]=$rev_str;
			$new_fields[$ABB_COL]=$fields[$ABB_COL] . "_$ref_abbr";
			$new_fields[$EQU_COL]=$new_rxn_eqn;
			$out_str=join "\t", @new_fields;
	
			print STDERR "\n";
			print STDERR "$tar_dir / $tar_rev\n";
			print STDERR "Original: $rxn_eqn\n";
			print STDERR "     New: $new_rxn_eqn\n";
			print STDERR "\n";

		}else{
			$out_str=$_;
		}	

		print OUT_FH "$out_str\n";
	}
}

close(IN_FH);
close(OUT_FH);

print STDERR "\n";
print STDERR "Done.\n";



























