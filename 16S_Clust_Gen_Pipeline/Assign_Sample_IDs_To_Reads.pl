#!/usr/bin/env perl

###############################################################################

use FindBin ();
use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_m $opt_o $opt_f $opt_r);

my $PIPELINE_UTIL_PATH="$FindBin::Bin/pipeline_utilities";
print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";
my $RENAME_BIN="$PIPELINE_UTIL_PATH/../../FASTA/Rename_FASTA_Records_with_MappingFile.pl";

getopts("m:o:fr");
my $usage = "usage: 

$0 

	-m <sample name to FASTA map file>
	-o <output read to groups/sample root filename>
	[-f (build multifasta file across groups flag)]
	[-r (overwrite existing files flag)]

	This script will generate a read to group/sample mapping file
	for Mothur, or whoever needs it.

	The input file is:
		<sample/group name> \\t <FASTA file> \\n
		...

	The output file is:
		<read id> \\t <sample/group name \\n
		...

	If the -r option is specified, the output files will overwrite exising
	files if they already exist.

	If the -f option is specified, then a multifasta file will be created
	which contains all the fasta sequences specified in the input sample to
	FASTA file.

	The output file(s) will be:
		<output filename root>.groups
		<output filename root>.fasta
		<output filename root>.read_mapping

	NOTE: Mothur doesn't like :, -, and \\'s in the read and sample ids, so this
	script will also convert all those in the read and sample ids to _'s.

	Also, to shorten read ID names, new read ID will be changed to <sample_id>_#####,
	where ##### is the index of that read.

";

if(!(
	defined($opt_m) && 
	defined($opt_o))){
	die $usage;
}

my $sample_to_fasta_fname=$opt_m;
my $output_fname=$opt_o;
my $build_mfasta=defined($opt_f)?1:0;
my $overwrite=defined($opt_r)?1:0;

print STDERR "\n";
print STDERR "Build multifasta? $build_mfasta\n";
print STDERR "Overwrite Existing? $overwrite\n";

###############################################################################
# Read in map file

print STDERR "\n";
print STDERR "Loading Sample to FASTA list: $sample_to_fasta_fname\n";
print STDERR "\n";

open(MAP_FH, "<$sample_to_fasta_fname") || die "Could not open $sample_to_fasta_fname\n";

my @sample_name_arr;
my @fasta_fname_arr;
my $num_maps_loaded=0;

while(<MAP_FH>){
	chomp;

	my ($sample_id, $fasta_fname)=split "\t";

	push @sample_name_arr, $sample_id;
	push @fasta_fname_arr, $fasta_fname;

	$num_maps_loaded++;
}

close(MAP_FH);

print STDERR "Number of map entries loaded: $num_maps_loaded\n";
print STDERR "\n";

###############################################################################

sub get_read_ids_from_fasta{
	my $fasta=shift;
	my $num_ids_in_fasta=0;

	open(FASTA_FH, "<$fasta") || die "Could not open $fasta.\n";

	my @ids;
	while(<FASTA_FH>){
		if($_=~/^>/){
			chomp;

			my $defline=$_;
			if($defline=~/^>(\S+)/){
				push @ids, $1;
				$num_ids_in_fasta++;
			}
		}
	}

	close(FASTA_FH);

	print STDERR "  Num IDs found in FASTA: $num_ids_in_fasta\n";
	return(\@ids);
}

#-------------------------------------------------------------------------------

sub output_sample_read_ids{
	my $filename=shift;
	my $group_name=shift;
	my $read_id_arr_ref=shift;
	my $num_lines_out=0;

	open(FH, ">>$filename") || die "Could not append to $filename.\n";

	foreach my $read_id (@{$read_id_arr_ref}){
		print FH "$read_id\t$group_name\n";
		$num_lines_out++;
	}

	close(FH);

	print STDERR "  Num Lines Output: $num_lines_out\n";

}

#-------------------------------------------------------------------------------

sub rename_id{
	my $read_id_arr_ref=shift;
	my $sample_id=shift;
	my $rename_map_hash_ref=shift;

	my @renamed_read_id_arr;

	my $i=0;
	foreach my $read_id (@{$read_id_arr_ref}){
		my $new_id=sprintf("$sample_id\_%06i", $i);
		push @renamed_read_id_arr, $new_id;
		${$rename_map_hash_ref}{$read_id}=$new_id;
		$i++;
	}

	return(\@renamed_read_id_arr);
}

###############################################################################

my $groups_fname="$output_fname.groups";
print STDERR "Output Groups Filename: $groups_fname\n\n";
if(-e $groups_fname){
	if(!$overwrite){
		die "Error: $groups_fname already exists.\n\n";
	}else{
		truncate $groups_fname, 0;
	}
}

# Generate read to group/sample mapping
my $total_reads_output=0;
my %mapping_hash;
for(my $i=0; $i<$num_maps_loaded; $i++){
	my $sample_id=$sample_name_arr[$i];
	my $fasta_fname=$fasta_fname_arr[$i];

	print STDERR "Cleaning Sample ID names. (Substituting :,-, and \\ with _'s.)\n";
	$sample_id=~s/:/_/g;
	$sample_id=~s/-/_/g;
	$sample_id=~s/\//_/g;
	
	print STDERR "Extracting Read IDs from $fasta_fname...\n";
	my $read_id_arr_ref=get_read_ids_from_fasta($fasta_fname);

	print STDERR "Renaming Read IDs...\n";
	my $renamed_read_id_arr_ref=rename_id($read_id_arr_ref, $sample_id, \%mapping_hash);

	print STDERR "Outputing Mapping Info for $sample_id...\n";
	output_sample_read_ids($groups_fname, $sample_id, $renamed_read_id_arr_ref);

	$total_reads_output+=($#{$read_id_arr_ref}+1);

	print STDERR "\n";
}

print STDERR "Total Reads Processed: $total_reads_output\n";

# Out mapping file
my $output_mapping_fname="$output_fname.read_mapping";
open(MAP_FH, ">$output_mapping_fname") || die "Could not open $output_mapping_fname for writing.\n";
foreach my $key (sort keys %mapping_hash){
	print MAP_FH "$key\t$mapping_hash{$key}\n";
}
close(MAP_FH);

# Concatenate all the fasta files together
if($build_mfasta){

	print STDERR "Building multifasta across all samples...\n";
	my $mfasta_fname="$output_fname.fasta";

	if(-e $mfasta_fname){
		if(!$overwrite){
			die "Error: $mfasta_fname already exists.\n\n";
		}else{
			truncate $mfasta_fname, 0;
		}
	}
	
	for(my $i=0; $i<$num_maps_loaded; $i++){
		`cat $fasta_fname_arr[$i] >> $mfasta_fname`;
	}

	print STDERR "Remapping FASTA file IDs...\n";
	my $ren_cmd_str="perl $RENAME_BIN -i $mfasta_fname -o $mfasta_fname.tmp -m $output_mapping_fname";
	print STDERR "'$ren_cmd_str'\n";
	print STDERR `$ren_cmd_str`;

	system("mv $mfasta_fname.tmp $mfasta_fname");
}

###############################################################################

print STDERR "done.\n";

