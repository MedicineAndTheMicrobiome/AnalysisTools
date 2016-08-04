#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_b $opt_p);
use FindBin;
getopts("f:b:p");

my $TAXA_LIST_EXT="taxa_list";
my $SUMMARY_TABLE_XLS_EXT="summary_table.xls";
my $R_INPUT_TXT_EXT="r_input.txt";
my $DEFAULT_BERGEY="$FindBin::Bin/bergey_ranks";
my $PAD_UNCLASSIFIED_BIN="$FindBin::Bin/Pad_Unclassified.pl";

my $usage = "usage:
$0
	-f <sample id to file name mapping file>
	[-p flag to pad unknown classifications]
	[-b <bergey classification file>]

	Sample ID to _countcum.txt Mapping File (-f):

	    <sample_id>\\t<*_countcum.txt file>\\n
	    <sample_id>\\t<*_countcum.txt file>\\n
	    <sample_id>\\t<*_countcum.txt file>\\n
		...
	
	Note: 
		1.) The *countcum.txt file should be a path of where this program can find each file.

	Taxonomic ranks automatically pulled:
		1.) Domain 
		2.) Phylum
		3.) Class 
		4.) Order 
		5.) Family
		6.) Genus

	Output Root (-o):
		The following output will be generated:
		<output root>.<taxonomic rank>.$TAXA_LIST_EXT         (list of all taxas used)
		<output root>.<taxonomic rank>.$SUMMARY_TABLE_XLS_EXT (tab separated text values for Excel)
		<output root>.<taxonomic rank>.$R_INPUT_TXT_EXT       (simplified space separated values for R)

	Default BergeyFile: $DEFAULT_BERGEY

";

if(!(
	defined($opt_f)
)){
	die $usage;
}

my $input_id_file_mapping_file=$opt_f;
my $bergey_class_file;
my $pad_unknown_classifications;

if(!defined($opt_b)){
	$bergey_class_file=$DEFAULT_BERGEY;
}else{
	$bergey_class_file=$opt_b;
}

if(!defined($opt_p)){
	$pad_unknown_classifications=0;
}else{
	$pad_unknown_classifications=1
}

my $output_file_root=$input_id_file_mapping_file;

if($pad_unknown_classifications){
	$output_file_root.=".padded";
}else{
	$output_file_root.=".unpadded";
}

my ($sample_arr_ref, $path_hash_ref)=load_samples_info($input_id_file_mapping_file);
my $taxon_hash_ref=load_bergey_class_file($bergey_class_file);
my %taxon_hash=%{$taxon_hash_ref};

#------------------------------------------------------------------------------

my %count_hash;
my %all_taxa_hash;
my %taxa_level_counts;

foreach my $sample_id(@{$sample_arr_ref}){
	my $path=${$path_hash_ref}{$sample_id};

	print STDERR "Working on $sample_id\n";

	my $counts_hash_ref=load_count_cumu($path, $pad_unknown_classifications);
	
	foreach my $taxa (keys %{$counts_hash_ref}){
		my $rank=${$taxon_hash_ref}{$taxa};
		my $count=${$counts_hash_ref}{$taxa};

		if(!defined($rank)){
			
			my @taxa_arr=split / /, $taxa;
			my $finest_grain=pop @taxa_arr;
			#print STDERR "$taxa: $finest_grain\n";
			if($finest_grain=~/\#uncertain_(.+)/){
				$rank=$1;
				my ($name, $phylum)=split /#/, $taxa;
				$taxa=$name;
			}else{
				print STDERR "Unknown rank : '$taxa'\n";
			}
		}
		#print STDERR "$rank: $taxa: ${$counts_hash_ref}{$taxa}\n";
		$count_hash{$sample_id}{$rank}{$taxa}=$count;
		$all_taxa_hash{$rank}{$taxa}=1;

		$taxa_level_counts{$rank}{$taxa}+=$count;
	}
}

my @ranks_of_interest=("domain", "phylum", "class", "order", "family", "genus");

foreach my $target_rank(@ranks_of_interest){

	my @taxa_arr=sort keys %{$all_taxa_hash{$target_rank}};
	my $outputfile;

	#------------------------------------------------------------------------------
	# Write out the taxa list

	open (TAXA_LIST_FH, ">$output_file_root\.$target_rank\.$TAXA_LIST_EXT") ||
		die "Could not open $outputfile\.$target_rank\.TAXA_LIST_EXT\n";
	foreach my $taxa(@taxa_arr){
		print TAXA_LIST_FH "$taxa\t$taxa_level_counts{$target_rank}{$taxa}\n";
	}
	close(TAXA_LIST_FH);

	#------------------------------------------------------------------------------
	# Write out the summary table

	open (SUMMARY_TABLE_XLS_FH, ">$output_file_root\.$target_rank\.$SUMMARY_TABLE_XLS_EXT") || 
		die "Could not open $outputfile\.$target_rank\.SUMMARY_TABLE_XLS_EXT\n";

	# Add up all the taxa hits
	my %sample_totals_hash;
	foreach my $taxa (@taxa_arr){
		foreach my $sample_id(@{$sample_arr_ref}){
			$sample_totals_hash{$sample_id}+=$count_hash{$sample_id}{$target_rank}{$taxa};
		}
	}

	# print out the column title/headers
	print SUMMARY_TABLE_XLS_FH (join "\t", ("Sample", "Total")) . "\t";
	print SUMMARY_TABLE_XLS_FH ((join "\t", @taxa_arr) . "\n");

	# print out the hits for each sample
	foreach my $sample_id(@{$sample_arr_ref}){
		my @val;
		foreach my $taxa(@taxa_arr){

			if(!defined($count_hash{$sample_id}{$target_rank}{$taxa})){
				push @val, 0;
			}else{
				push @val, $count_hash{$sample_id}{$target_rank}{$taxa};
			}
		}
		my $out_str=join "\t", ($sample_id, $sample_totals_hash{$sample_id}, @val);
		print SUMMARY_TABLE_XLS_FH "$out_str\n";
	}

	close(SUMMARY_TABLE_XLS_FH);

	#------------------------------------------------------------------------------
	# Output R

	open (R_INPUT_TXT_FH, ">$output_file_root\.$target_rank\.$R_INPUT_TXT_EXT") || die "Could not open $outputfile\.$target_rank\.R_INPUT_TXT_EXT\n";

	# Dump column titles
	my @pseudo_titles;
	push @pseudo_titles, "";
	foreach(my $i=0; $i<=$#taxa_arr; $i++){
		my @taxa_parts=split / /, $taxa_arr[$i];
		push @pseudo_titles, pop(@taxa_parts);
	}
	my $outstr=join " ", @pseudo_titles;
	print R_INPUT_TXT_FH "$outstr\n";

	# print out the hits for each sample
	foreach my $sample_id(@{$sample_arr_ref}){
		my @val;
		foreach my $taxa(@taxa_arr){
			if(!defined($count_hash{$sample_id}{$target_rank}{$taxa})){
				push @val, 0;
			}else{
				push @val, $count_hash{$sample_id}{$target_rank}{$taxa};
			}
		}

		# Change all spaces in the sample id to underscores
		$sample_id=~s/ /_/g;

		my $out_str=join " ", ($sample_id,  @val);
		print R_INPUT_TXT_FH "$out_str\n";
	}

	close(R_INPUT_TXT_FH);

	#------------------------------------------------------------------------------
}

print STDERR "Done.\n\n";

###############################################################################

sub load_count_cumu{
	my $path=shift;
	my $pad_unclassified_counts=shift;

	my %count_hash;

	if($pad_unclassified_counts){
		open(COUNTS, "$PAD_UNCLASSIFIED_BIN < $path |") || die "Could not open $path for reading.\n";
	}else{
		open(COUNTS, "$path") || die "Could not open $path for reading.\n";
	}

	while(<COUNTS>){
		chomp;

		#print "before: $_\n";

		$_=~s/"//g;
		my @line=split /\t/, $_;

		my $total=shift @line;
		my $leaf=pop @line;
		if($leaf ne "(any)"){
			push @line, $leaf;
		}

		my $taxa=join " ", @line;
		$count_hash{$taxa}=$total;

		#print "after: $total / $taxa\n\n";

	}

	close(COUNTS);

	return(\%count_hash);
}

###############################################################################

sub load_bergey_class_file{
	my $bergey_class_file=shift;

	my %class_hash;

	open(FH, "<$bergey_class_file") || die "Could not open $bergey_class_file for reading.\n";

	while(<FH>){
		chomp;
		my ($rank, $name)=split /\t/, $_;
		$name=~s/"//g;
		$class_hash{$name}=$rank;
	}

	close(FH);

	# Stick in this pseudo record in case the classification is unknown
	$class_hash{"Unknown"}="domain";

	return(\%class_hash);
}


###############################################################################

sub load_samples_info{
	my $input_id_file_mapping_file=shift;
	my @sample_list;
	my %path_hash;

	open(FH, "<$input_id_file_mapping_file") || die "Could not open sample mapping list: $input_id_file_mapping_file\n";

	while(<FH>){
		chomp;
		my ($sample_id, $file_path)=split /\t/, $_;

		if(!defined($sample_id) || $sample_id eq ""){
			die "Error loading sample mapping list, null Sample ID: \"$_\"\n";	
		}elsif(!defined($file_path) || $file_path eq ""){
			die "Error loading sample mapping list, null File Path: \"$_\"\n";	
		}

		push @sample_list, $sample_id;
		$path_hash{$sample_id}=$file_path;
	}
	
	return(\@sample_list, \%path_hash);
}


