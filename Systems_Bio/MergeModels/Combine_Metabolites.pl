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
		-a <Input metabolites, model A>
		-b <Input metabolites, model B>	
		-o <Output filename root>

	This script will read in two metabolite files and merge them,
	removing duplicates. Identifying duplicates are based on the
	metabolite ID.  The description from Model A is taken over 
	model B.  Compartments are also merged.
	
";

if(!defined($opt_o) || 
	!defined($opt_a) || !defined($opt_b)
){
	die $usage;
}

my $MetFileA=$opt_a;
my $MetFileB=$opt_b;
my $Output=$opt_o;

###############################################################################

print STDERR "Met File A: $MetFileA\n";
print STDERR "Met File B: $MetFileB\n";
print STDERR "Output: $Output\n";
print STDERR "\n";

###############################################################################

my $ID_POS=0;
my $DESC_POS=1;
my $COMP_POS=2;

###############################################################################

sub load_metabolites{
	my $fname=shift;
	my $hash_ref=shift;
	
	print STDERR "Loading: $fname\n";
	open(IN_FH, "<$fname") || die "Could not open $fname for reading.\n";

	while(<IN_FH>){
		chomp;
		
		#print("$_\n");
		my ($id, $desc, $comp)=split "\t", $_;

		# Skip the header line
		if($id eq "abbreviation"){
			next;
		}

		# Clean up compartment string
		$comp=~s/"//g;
		$comp=~s/ //g;

		# Get old record
		my $prior_met_rec=${$hash_ref}{$id};

		# Replace old record with new record, but append compartment info
		if(defined($prior_met_rec)){
			my ($prior_desc, $prior_comp)=split "\t", $prior_met_rec;

			#print STDERR "'$prior_comp'/'$comp'\n";
			
			if($prior_comp eq "" && $comp eq ""){
				$comp="";
			}elsif($prior_comp ne "" && $comp eq ""){
				$comp=$prior_comp;
			}elsif($prior_comp eq "" && $comp ne ""){
				$comp=$comp;
			}elsif($prior_comp ne "" && $comp ne ""){
				$comp="$prior_comp,$comp";	
			}

			#print STDERR "$comp\n\n";
				
		}

		# Save new record
		${$hash_ref}{$id}="$desc\t$comp";
	}

	close(IN_FH);

}

###############################################################################

sub dump_metabolites{
	my $fname=shift;
	my $hash_ref=shift;

	print STDERR "Outputing: $fname\n";
	open(OUT_FH, ">$fname") || die "Could not open $fname for writing.\n";
	
	print OUT_FH "abbreviation\tname\tcompartment\n";	

	foreach my $id(sort keys %{$hash_ref}){
		my ($desc, $comps)=split "\t", ${$hash_ref}{$id};
		
		# Make unique compartment list
		my %hash;
		my @comp_arr=split ",", $comps;
		foreach my $comp_unit(@comp_arr){
			$hash{$comp_unit}=1;
		}
		@comp_arr=sort keys %hash;
		$comps=join ",", @comp_arr;

		# Output metabolite record
		print OUT_FH "$id\t$desc\t$comps\n";
	}

	close(OUT_FH);
}

###############################################################################

my %met_hash;

load_metabolites($MetFileB, \%met_hash);
load_metabolites($MetFileA, \%met_hash);
dump_metabolites($Output, \%met_hash);

###############################################################################

print STDERR "\n";
print STDERR "Done.\n";

















