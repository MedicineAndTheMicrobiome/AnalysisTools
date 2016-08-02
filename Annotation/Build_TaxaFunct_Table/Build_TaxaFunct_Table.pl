#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_t $opt_f $opt_n $opt_c $opt_x $opt_o);

getopts("t:f:n:c:x:o:");

my $usage = "
usage:

$0
	-t <taxa list>
	-f <function list>

	-n <sample annotation tsv file>
	-c <function column, from columns 4-9>
	-x <sample taxa_ids tsv file	

	-o <output filename root>

	This script will read in the sample annotation and 
	sample taxa ids file and generate a contigency table
	based on the specified taxa and function.

	
	The sample taxa_ids tsv file (-x) should contain the following
	columns.
		1.) Read id
		2.) Taxa assigned to read
		3.) Highest assigned taxa id (i.e. Superkingdom)
		4.) Next high assigned taxa id (e.g. Phylum)
		...
		10.) Lowest asigned taxa id (e.g. Genus)

	The sample annotation tsv file (-n) should contain the following
	columns.  This file should be based on the accumulated evidence
	across all the hits.  Columns 4-9 should be specified inside the
	file specified with -f.
		1.) Read ID
		2.) Percent Composite Identity
		3.) Length
		4.) Description
		5.) PFam IDs
		6.) TIGRFam IDs
		7.) GO Process
		8.) GO Function
		9.) EC

";

if(
	!defined($opt_t) || 
	!defined($opt_f) || 
	!defined($opt_n) || 
	!defined($opt_c) || 
	!defined($opt_x) || 
	!defined($opt_o) 
){
	die $usage;
}

###############################################################################

my $taxa_list_fn=$opt_t;
my $function_list_fn=$opt_f;
my $samp_anno_fn=$opt_n;
my $samp_anno_col=$opt_c;
my $samp_taxa_fn=$opt_x;
my $output_fn_root=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "Target Taxa Filename: $taxa_list_fn\n";
print STDERR "Target Function Filename: $function_list_fn\n";
print STDERR "Sample Annotation Filename: $samp_anno_fn\n";
print STDERR "\tTargeted Column: $samp_anno_col\n";
print STDERR "Sample Taxa Filename: $samp_taxa_fn\n";
print STDERR "Output Filename Root: $output_fn_root\n";
print STDERR "\n";

###############################################################################

sub load_list{
	my $fname=shift;
	my @list;
	my %dup;

	print STDERR "Loading list: $fname\n";
	open(FH, "<$fname") || die "Could not open list: $fname\n";
	while(<FH>){
		chomp;
		my ($item, $junk)=split "\t", $_;

		if(!defined($dup{$item})){
			push @list, $item;
		}else{
			print STDERR "Warning: $item was removed.  It was duplicated.\n";
		}

		$dup{$item}=1;
	}

	return(\@list);
}

sub print_list{
	my $arr_ref=shift;
	for(my $i=0; $i<=$#{$arr_ref}; $i++){
		print STDERR "${$arr_ref}[$i]\n";
	}
}

###############################################################################

sub load_function{
	my $func_fn=shift;
	my $tar_fun_arr_ref=shift;
	my $col=shift;

	$col--; # Count from 0

	# Load targets into hash
	my %target_fun_hash;
	foreach my $fun_id (@{$tar_fun_arr_ref}){
		$target_fun_hash{$fun_id}=1;
	}

	# Keep hash
	my %keep_hash;
	my $target_ct=0;
	my $undefined_fun_ct=0;
	my $other_fun_ct=0;

	open(FH, "<$func_fn") || die "Could not open $func_fn\n";

	while(<FH>){
		chomp;
		my @columns=split "\t", $_;
		my $read_id=$columns[0];
		my $tar_fun=$columns[$col];

		if($tar_fun eq "" || !defined($tar_fun)){
			$undefined_fun_ct++;
			$keep_hash{$read_id}="UndefinedFunction";
		}else{
			if(defined($target_fun_hash{$tar_fun})){
				$keep_hash{$read_id}=$tar_fun;
				$target_ct++;
			}else{
				$other_fun_ct++;
				$keep_hash{$read_id}="OtherFunction";
			}
		}
	}

	close(FH);

	print STDERR "Functions Identified:\n";
	print STDERR "  Targeted Found: $target_ct\n";
	print STDERR "  Undefined Function: $undefined_fun_ct\n";
	print STDERR "  Other Function: $other_fun_ct\n";

	return(\%keep_hash, $target_ct, $undefined_fun_ct, $other_fun_ct);
}

sub load_taxa{
	my $taxa_fn=shift;
	my $tar_taxa_arr_ref=shift;
	
	# Load targets into hash
	my %target_tax_hash;
	foreach my $tax_id (@{$tar_taxa_arr_ref}){
		$target_tax_hash{$tax_id}=1;
	}

	# Keep hash
	my %keep_hash;
	my $target_ct=0;
	my $undefined_tax_ct=0;
	my $other_tax_ct=0;
	
	open(FH, "<$taxa_fn") || die "Could not open $taxa_fn\n";

        while(<FH>){
                chomp;
                my @columns=split "\t", $_;

		my $read_id=$columns[0];
		my $taxa_hit=0;
		my $not_nulls=0;

		for(my $i=1; $i<=$#columns; $i++){
			my $tax_id=$columns[$i];
			if($tax_id eq "NULL"){

			}else{
				$not_nulls++;

				if(defined($target_tax_hash{$tax_id})){
					$keep_hash{$read_id}=$tax_id;
					$taxa_hit++;
				}
			}
		}

		if($taxa_hit>0){
			$target_ct++;
		}elsif($not_nulls>0){
			$other_tax_ct++;
			$keep_hash{$read_id}="OtherTaxa";
		}else{
			$undefined_tax_ct++;
			$keep_hash{$read_id}="UnknownTaxa";
		}

	}

        close(FH);

        print STDERR "Taxa Identified:\n";
        print STDERR "  Targeted Found: $target_ct\n";
        print STDERR "  Undefined Taxa: $undefined_tax_ct\n";
        print STDERR "  Other Taxa: $other_tax_ct\n";

        return(\%keep_hash, $target_ct, $undefined_tax_ct, $other_tax_ct);

}

sub build_contigency_table{
	my $func_hash_ref=shift;
	my $taxa_hash_ref=shift;
	my $func_list_ref=shift;
	my $taxa_list_ref=shift;

	my %ids;

	foreach my $id(keys %{$func_hash_ref}){
		$ids{$id}=1;
	}

	foreach my $id(keys %{$taxa_hash_ref}){
		$ids{$id}=1;
	}

	my %contigency_table;

	# Initialize
	my @func_list=@{$func_list_ref};
	my @taxa_list=@{$taxa_list_ref};
	push @func_list, ("OtherFunction", "UndefinedFunction");
	push @taxa_list, ("OtherTaxa", "UnknownTaxa");
	foreach my $func(@func_list){
		foreach my $taxa(@taxa_list){
			${$contigency_table{$taxa}}{$func}=0;
		}
	}

	# Perform counts
	my $tot=0;
	foreach my $id(keys %ids){

		my $func=${$func_hash_ref}{$id};
		my $taxa=${$taxa_hash_ref}{$id};

		if(!defined($func)){
			$func="OtherFunction";
		}

		if(!defined($taxa)){
			$taxa="OtherTaxa";
		}


		my $before=${$contigency_table{$taxa}}{$func};
		my $after=$before+1;

		${$contigency_table{$taxa}}{$func}=$after;
		$tot++;

	}
	print STDERR "Total Reads Analyzed = $tot\n";

	return(\%contigency_table);
}


sub print_contigency_table{
	my $fname=shift;
	my $func_list_ref=shift;
	my $taxa_list_ref=shift;
	my $contigency_table_ref=shift;

	my @func_list=@{$func_list_ref};
	my @taxa_list=@{$taxa_list_ref};

	push @func_list, ("OtherFunction", "UndefinedFunction");
	push @taxa_list, ("OtherTaxa", "UnknownTaxa");

	open(FH, ">$fname") || die "Could not open $fname.\n";

	print FH "TaxaIDs";
	foreach my $func(@func_list){
		print FH "\t$func";
	}
	print FH "\n";

	my $tot=0;
        foreach my $taxa(@taxa_list){
		print FH "$taxa";
		foreach my $func(@func_list){
			print FH "\t";
                        print FH ${$contigency_table_ref}{$taxa}{$func};
			$tot+=${$contigency_table_ref}{$taxa}{$func};
                }
		print FH "\n";
        }
	close(FH);

	print STDERR "Total Counts In Contigency Table = $tot\n";
}


##############################################################################

my $taxa_list_ref=load_list($taxa_list_fn);
my $func_list_ref=load_list($function_list_fn);

print STDERR "\n";
print STDERR "Function Targets:\n";
print_list($func_list_ref);
print STDERR "\n";
print STDERR "Taxa Targets:\n";
print_list($taxa_list_ref);
print STDERR "\n";

###############################################################################

my ($kept_reads_funct_hash, $num_kept, $num_undef, $num_other)=
	load_function($samp_anno_fn, $func_list_ref, $samp_anno_col);

print STDERR "\n";

my ($kept_reads_taxa_hash, $num_kept, $num_undef, $num_other)=
	load_taxa($samp_taxa_fn, $taxa_list_ref);

print STDERR "\n";

my $contig_table_ref=build_contigency_table(
	$kept_reads_funct_hash, $kept_reads_taxa_hash,
	$func_list_ref, $taxa_list_ref
);
	
print_contigency_table(
	"$output_fn_root\.ctgcy.tsv", 
	$func_list_ref,
	$taxa_list_ref,
	$contig_table_ref
);	


###############################################################################

print STDERR "Done.\n";
