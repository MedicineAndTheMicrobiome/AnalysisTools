#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_p $opt_a $opt_t $opt_o $opt_f $opt_c $opt_e);

getopts("p:a:t:o:f:c:e");

my $usage = "
usage:
    $0

	Either:
	 -p <annotation/taxa/output 3-tuple assignment list file>

	or

	All 3:
	 -a <annotation assignment file>
	 -t <taxa assignment file>
	 -o <output filtered annotation assignment root>

	Required:
	 -f <target function list>
	 -c <target annotation column, read id is column 0 (first)>

	Optional:
	[-e (exact match the target function list)]

	All files should have their fields tab separated.

	Assignment list file:
		annotation filename, taxonomy filename,  output filename

	Annotation file (-a):
		query_id, ...<remaining annotation>
		
		where the remaining annotation columns are:
			3: Description
			4: PFam
			5: TIGRFam
			6: GO Process 
			7: GO Function
			8: EC

	Taxa assignment file (-t):
		query_id, taxa_id,
                superkingdom, phylum, class, order, family, genus, species

	Target list of function ID to exact match (-f)
		function_1
		function_2
		...
		function_N

	When the -e option is specified an exact match will be used
	for the annotation function list.  By default, if the annotation
	for a read is a semi-colon separated list of functions, the 
	list will be split into its items and if any of the items match
	the target, the read will be selected.  With the -e option,
	the list will not be split and must match exactly with the
	targeted function.

	This script will read in the target function list and inside the
	specified column, extract all the read IDs that match.  Then it
	will extract out all the reads in the taxa assigment file.  The
	result is a new taxa assignment file which can be use to create
	a summary table.

";

###############################################################################

if(
	!defined($opt_f) ||
	!defined($opt_c)
){
	die $usage;
}

if(!defined($opt_p) && (!defined($opt_a) || !defined($opt_t))){
	die "You need to specify either the -p option or the -a and -t option\n";
}

my $AnnotationTaxaPairList=$opt_p;
my $AnnotationFile=$opt_a;
my $TaxaFile=$opt_t;
my $TargetFunctionListFile=$opt_f;
my $TargetColumn=$opt_c;
my $OutputFilenameRoot=$opt_o;
my $ExactMatch=defined($opt_e)? 1:0;

my $use_list;

if(defined($opt_p)){
	$use_list=1;		
}elsif(defined($opt_a) && defined($opt_t) && defined($opt_o)){
	$use_list=0;
}else{
	print STDERR "Need to specify a list of the a/t/o option.\n";
}

if($use_list){
	print STDERR "Annotation/Taxa Pair List: $AnnotationTaxaPairList\n";
}else{
	print STDERR "Annotation File: $AnnotationFile\n";
	print STDERR "Taxonomy File: $TaxaFile\n";
	print STDERR "Output Filename Root: $OutputFilenameRoot\n";
}
print STDERR "\n";
print STDERR "Target Function List: $TargetFunctionListFile\n";
print STDERR "Target Column: $TargetColumn\n";
print STDERR "Perform exact match? $ExactMatch\n";

###############################################################################

sub load_hash{

	my $fname=shift;
	my %hash;
	
	open(FH, "<$fname") || die "Could not open $fname\n";

	while(<FH>){
		chomp;
		my ($item, $comment)=split "\t", $_;
		$hash{$item}=1;
	}
	
	close(FH);

	if(defined($hash{""})){
		delete $hash{""};
	}

	my $num_functions=keys %hash;
	print STDERR "Num Functions to Target: $num_functions\n";

	return(\%hash);	
}

#------------------------------------------------------------------------------

sub get_ids_from_annotation{

	my $filename=shift;
	my $column=shift;
	my $target_function_hash=shift;
	my $exact_match=shift;

	open(ANNOFH, "<$filename") || die "Could not open Annotation File: $filename\n";

	my %read_id_hash;
	while(<ANNOFH>){
		chomp;
		my @columns=split "\t", $_;

		my $cur_ann=$columns[$column];
		my $cur_id=$columns[0];
	
		my @items;
		if($exact_match){
			push @items, $cur_ann;
		}else{
			@items=split ";", $cur_ann;
		}
		
		foreach my $item(@items){
			if(defined(${$target_function_hash}{$item})){
				$read_id_hash{$cur_id}=1;	
			}
		}

		if(defined(${$target_function_hash}{$cur_ann})){
			$read_id_hash{$cur_id}=1;	
		}

	}

	close(ANNOFH);

	my $num_read_ids=keys %read_id_hash;
	print STDERR "Num Read IDs found with targeted annotation: $num_read_ids\n";
	
	return(\%read_id_hash);
}

#------------------------------------------------------------------------------

sub extract_taxa_by_read_id{

	my $filename=shift;
	my $output_filename=shift;
	my $target_read_ids_hash_ref=shift;

	open(FH, "<$filename") || die "Could not open Taxa File: $filename\n";
	open(OUT_FH, ">$output_filename") || die "Could not open Output Taxa File: $output_filename\n";

	while(<FH>){
		chomp;

		my @columns=split "\t", $_;
		if(defined(${$target_read_ids_hash_ref}{$columns[0]})){
			print OUT_FH "$_\n";
		}	
	}	
	
	close(FH);
	close(OUT_FH);

}

###############################################################################
# Load functions to keep read ids for

my $function_filt_hash_ref=load_hash($TargetFunctionListFile);

print STDERR "Targeted Functions:\n";
foreach my $funct (keys %{$function_filt_hash_ref}){
	print STDERR "'$funct'\n";
}
print STDERR "\n";

###############################################################################

if($use_list){
	open(ANFH, "<$AnnotationTaxaPairList") || die "Could not open $AnnotationTaxaPairList\n";
	
	while(<ANFH>){
		chomp;
		my ($function_fn, $taxa_fn, $output_fn)=split "\t", $_;

		print STDERR "\n";
		print STDERR "Working on:\n";
		print STDERR "   Funct: $function_fn\n";
		print STDERR "   Taxa: $taxa_fn\n";
		print STDERR "   Output: $output_fn\n";

		if($function_fn eq "" || $taxa_fn eq "" || $output_fn eq ""){
			die "Error: Column missing in annotation/taxa list.\n";
		}

		print STDERR "Getting read IDs matching target function...\n";
		my $read_id_hash_ref=get_ids_from_annotation($function_fn, $TargetColumn, $function_filt_hash_ref, $ExactMatch);

		print STDERR "Getting taxonomy based on reads IDs...\n";
		extract_taxa_by_read_id($taxa_fn, $output_fn, $read_id_hash_ref);
	
		print STDERR "ok.\n\n";
	}

	close(ANFH);

}else{
	my $read_id_hash_ref=get_ids_from_annotation($AnnotationFile, $TargetColumn, $function_filt_hash_ref, $ExactMatch);
	extract_taxa_by_read_id($TaxaFile, $OutputFilenameRoot, $read_id_hash_ref);
}

###############################################################################

print STDERR "Done.\n";

