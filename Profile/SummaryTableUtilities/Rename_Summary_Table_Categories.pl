#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_m $opt_o $opt_r $opt_l);
getopts("i:m:o:rl:");

my $usage = "usage:
$0
		-i <input summary_table.tsv name>
		-m <mapping file>
		[-o <output summary_table.tsv name>]
		[-l \"<list separator>\"]
		[-r (replace old name, instead of placing it in parenthesis.)]
	
	Translates/changes the first row of the summary table (the taxa/categories names)
	from to a new name.

	The mapping file should be:
		<old name>\\t<new name>\\n

	If the -l option is used, the string will be split into a list, translated,
	and then reassembled.

";

if(!(
	defined($opt_i) || 
	defined($opt_m)
)){
	die $usage;
}

my $input_file_name=$opt_i;
my $map_name=$opt_m;
my $replace;

if(defined($opt_r)){
	$replace=1;
}else{
	$replace=0;
}

my $list_sep;
if(defined($opt_l)){
	$list_sep=$opt_l;
}else{
	$list_sep=undef;
}

my $output_file_name;
if(defined($opt_o)){
    	$output_file_name=$opt_o;
}else{
	$output_file_name=$input_file_name;
	$output_file_name=~s/\.summary_table\.tsv$//;
	$output_file_name=~s/\.summary_table\.xls$//;
	$output_file_name="$output_file_name\.mapped.summary_table.tsv";
}

print STDERR "Input Filename: $input_file_name\n";
print STDERR "Output Filename: $output_file_name\n";
print STDERR "GO Mapping: $map_name\n";
print STDERR "Replace?: $replace\n";

print STDERR "List separator: \"$list_sep\"\n";
		
##################################################################################

sub load_mappings{
	my $map_fname=shift;
	my %map_hash;

	open(FH, "<$map_fname") || die "Could not open $map_fname\n";
	while(<FH>){
		chomp;
		my ($id, $new)=split "\t", $_;
		$map_hash{$id}=$new;
	}
	close(FH);
	return(\%map_hash);
}

##################################################################################

open(OUT_FH, ">$output_file_name") || die "Could not open > $output_file_name";

my $map_hash_ref=load_mappings($map_name);

open(ST_FH, "<$input_file_name") || die "Could not open < $input_file_name\n";

my $header=<ST_FH>;
chomp $header;

my @header_fields=split "\t", $header;

my $sample_id_field=shift @header_fields;
my $counts_field=shift @header_fields;

my @new_field_name;


sub translate{
	my $old=shift;
	my $replace=shift;
	
	my $new=${$map_hash_ref}{$old};
	if(!defined($new)){
		return($old);
	}else{
		if($replace){
			return "$new";
		}else{
			return "$new ($old)";
		}
	}
}

if(!defined($list_sep)){

	foreach my $old_name (@header_fields){
		push @new_field_name, translate($old_name, $replace);
	}

}else{

	foreach my $old_name (@header_fields){
	
		my @components=split $list_sep, $old_name;
		my @new_components;
	
		foreach my $old_comp(@components){
			push @new_components, translate($old_comp, $replace);	
		}		

		push @new_field_name, (join $list_sep, @new_components);
	
	}
	
}

my $new_header=join "\t", ($sample_id_field, $counts_field, @new_field_name);
print OUT_FH "$new_header\n";

##################################################################################
# Output remaining

while(<ST_FH>){
	print OUT_FH "$_";
}

print STDERR "Done.\n";

##################################################################################



