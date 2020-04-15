#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_c $opt_m $opt_l $opt_o $opt_r $opt_p $opt_h $opt_d $opt_t $opt_a $opt_C $opt_f);

getopts("i:c:m:orphdtaC:fl:");

my $usage = "
	usage:
	$0
		-i <Input file name>
		-m <Map file, tab-separated>
		
		[-l \"<within field list separator>\">]
	
		[-o (overwrite original file)]	
		[-c <list of Columns to translate, starting from 0, default all>]

		Action options
		[-r (Replace column, default)]
		[-p (aPpended column)]
		[-h (adHere translation next to original value)]

		Input file options
		[-d (auto Detected, default)]
		[-t (Tab delimited)]
		[-a (commA delimited)]

		Maps File Options
		[-C <old key>,<new value>]

		Optimization options
		[-f (Build hash out of ID's in specified column, then only load those IDs from the Map file)]

	Reads map file into hash, then for each line in the input file, the columns
	specified will be translated to the value in the map file.

	If a mapping can not be made for the column, then it is left alone.

	Map file should have the format, unless the -C option is specified:
		<key>\\t<new value>\\n

	If the -C option is specified, then then alternative columns will be selected for
	the key to new value.

	For example -C 1,0 will swap the key to the second column (1) and the first column (0) will be the
	new value.

	Output goes to STDOUT.

";

if(!defined($opt_i) || !defined($opt_m)){
	die $usage;
}

my $InputFile=$opt_i;
my @Columns=split /,/, $opt_c;
my $MapFile=$opt_m;

my $Overwrite;
if(defined($opt_o)){
        $Overwrite=1;
}else{
	$Overwrite=0;
}	

# Replace Action
my $Action="Rep";
if(defined($opt_r)){
	$Action="Rep"; # Replace
}
if(defined($opt_p)){
	$Action="App"; # Append
}
if(defined($opt_h)){
	$Action="Adh"; # Adhere
}


# Delimitors
my $Delimitor="";
if(defined($opt_d)){
	$Delimitor="";
}
if(defined($opt_t)){
	$Delimitor="\t";
}
if(defined($opt_a)){
	$Delimitor=",";
}

# Optimization
my $Preread=0;
if(defined($opt_f)){
	$Preread=1;
}

# Within field list separator
my $ListSep="";
if(defined($opt_l)){
	$ListSep=$opt_l;
}

my $AltColKey=0;
my $AltColNew=1;
if(defined($opt_C)){
	($AltColKey, $AltColNew)=split /,/, $opt_C;
}

###############################################################################

print STDERR "Input: $InputFile\n";
print STDERR "Map: $MapFile\n";
print STDERR "Action: $Action\n";
print STDERR "Delimitor: '$Delimitor'\n";
print STDERR "Overwrite Original: $Overwrite\n";
print STDERR "List Separator: $ListSep\n";

if($AltColKey!=-1){
	print STDERR "Alternative Key/Map Columns Specified: Key=$AltColKey Map=$AltColNew\n";
}

print STDERR "Columns to Map:\n";
if($#Columns==-1){
	print STDERR "  ALL\n";
}else{
	foreach my $col(@Columns){
		print STDERR "  $col\n";
	}
}

###############################################################################

if($Delimitor eq ""){
	my $num_commas=0;
	my $num_tabs=0;

	print STDERR "Attempting to autodetected delimitor.\n";
	open(IN_FH, "head -n 20 $InputFile|") || die "Could not open $InputFile.\n";
	while(<IN_FH>){
		chomp;
		my $test=$_;
		$num_commas+=($test=~tr/,/,/);
		$num_tabs+=($test=~tr/\t/\t/);
	}
	close(IN_FH);

	print STDERR "Num Tabs: $num_tabs\n";
	print STDERR "Num Commas: $num_commas\n";
	if($num_tabs>(1.2*$num_commas)){
		$Delimitor="\t";
	}elsif($num_commas>(1.2*$num_tabs)){
		$Delimitor=",";
	}else{
		print STDERR "Could not detect delimitor.\n";
		print STDERR "Maybe it doesn't matter. Assuming <tab>\n";
		$Delimitor="\t";
		
	}
	print STDERR "I believe the delimitor is '$Delimitor'\n";
	print STDERR "\n";
}

###############################################################################

my %map;

if($Preread){
	print STDERR "Prereading input file for the mappings to retain.\n";
	open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";
	while(<IN_FH>){
		chomp;
		my @array=split /$Delimitor/, $_;

		if($#Columns==-1){
			for(my $i=0; $i<=$#array; $i++){
				$map{$array[$i]}="";
			}
		}else{
			foreach my $col(@Columns){
				$map{$array[$col]}="";
			}
		}
	}
	close(IN_FH);
	print STDERR "Done.\n";

	#foreach my $item(keys %map){
	#	print STDERR "'$item'\n";
	#}
}

###############################################################################

my $num_map_entries=0;

print STDERR "Loading map file.\n";
open(MAP_FH, "<$MapFile") || die "Could not open $MapFile\n";

if(!$Preread){
	print STDERR "Loading entire map file into memory.\n"; 
	while(<MAP_FH>){
		chomp;
		my @fields=split /\t/, $_;

		my $key=$fields[$AltColKey];
		my $val=$fields[$AltColNew];

		if($val ne ""){
			$map{$key}=$val;
		}else{
			$map{$key}="UNMAPPED";
		}
		$num_map_entries++;
		if(!($num_map_entries % 1000000)){
			print STDERR ".";
		}
	}
}else{
	print STDERR "Loading subset of map file into memory.\n"; 
	while(<MAP_FH>){
		chomp;
		my @fields=split /\t/, $_;

		my $key=$fields[$AltColKey];
		my $val=$fields[$AltColNew];

		if(defined($map{$key})){
			if($val ne ""){
               			$map{$key}=$val;
                	}else{
                        	$map{$key}="UNMAPPED";
                	}
			$num_map_entries++;
				if(!($num_map_entries % 1000)){
				print STDERR ".";
			}
		}
	}

}
close(MAP_FH);

print STDERR "\nNum Map Entries Read: $num_map_entries\n";

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;
open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

sub remap{
	my $action=shift;
	my $value=shift;
	my $delim=shift;
	my $list_separator=shift;
	
	my $remap="";

	my @terms=();
	if($list_separator ne ""){
		@terms=split /$list_separator/, $value;
	}else{
		@terms=($value);
	}

	my @remapped_terms;
	foreach my $term(@terms){
		if(defined($map{$term})){
			if($Action eq "Rep"){
				$remap=$map{$term};
			}elsif($Action eq "App"){
				$remap=$term . $delim . $map{$term};
			}elsif($Action eq "Adh"){
				$remap=$term . " (" . $map{$term} . ")";
			}else{
				print STDERR "Unknown action: $Action\n";
			}
		}else{
			$remap=$term;
		}
		push @remapped_terms, $remap;
	}

	my $remapped_str=join $list_separator, @remapped_terms;

	return($remapped_str);
}


my $fh;
my $temp_filename;
if($Overwrite){
	($fh, $temp_filename)=File::Temp::tempfile();
	print STDERR "Temp Filename: $temp_filename\n";
}

while(<IN_FH>){
	chomp;
	my @array=split /$Delimitor/, $_;

	if($#Columns==-1){
		for(my $i=0; $i<=$#array; $i++){
			$array[$i]=remap($Action, $array[$i], $Delimitor, $ListSep);
		}
	}else{
		foreach my $col(@Columns){
			$array[$col]=remap($Action, $array[$col], $Delimitor, $ListSep);
		}
	}

	my $outstr=join "$Delimitor", @array;

	if($Overwrite){
		print $fh "$outstr\n";
	}else{
		print "$outstr\n";
	}
	$input_lines=$input_lines+1;
}

close(IN_FH);

print STDERR "Num Lines Processed: $input_lines\n";

if($Overwrite){
	print STDERR "Overwriting original file.\n";
	`mv $temp_filename $InputFile`;
}

print STDERR "Done.\n";
