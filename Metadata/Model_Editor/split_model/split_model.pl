#!/usr/bin/env perl

###############################################################################

use strict;
use File::Basename;
use Getopt::Std;
use vars qw ($opt_m $opt_d);

getopts("m:d:");

my $usage = "
	usage:
	$0
		Input File Parameters:
		-m <input model, .mod file >
		-d <output destination, i.e. directory>

	This script will read in the model file and create a directory
	where the model's variable lists will be generated.  

";

if(
	!defined($opt_m) ||
	!defined($opt_d) 
){
	die $usage;
}

my $InputModelFile=$opt_m;
my $OutputDirectory=$opt_d;

print STDERR "Input Model File: $InputModelFile\n";
print STDERR "Output Directory: $OutputDirectory\n";
print STDERR "\n";

my ($model_root)=fileparse($InputModelFile);
$model_root=~s/\.mod$//;
print STDERR "Model Root: $model_root\n";
print STDERR "\n";


my $model_dir="$OutputDirectory/$model_root\.dir";
my $grp_dir="$OutputDirectory/$model_root\.dir/groups";

print "Output Model Dir: '$model_dir'\n";
print "  Groups Dir: '$grp_dir'\n";

mkdir $model_dir;
mkdir $grp_dir;

###############################################################################

sub write_list_to_file{
	my $arr_ref=shift;
	my $fname=shift;

	open(OUTPUT_FH, ">$fname") || die "Could not open output file: $fname\n";

	for(my $i=0; $i<=$#{$arr_ref}; $i++){
		print OUTPUT_FH "${$arr_ref}[$i]\n";
	}

	close(OUTPUT_FH);
}

###############################################################################
# Load the model

my %map_hash;

print STDERR "\n";
print STDERR "Loading Model File...\n";
open(MOD_FH, "<$InputModelFile") || die "Could not open map file: $InputModelFile\n";

while(<MOD_FH>){
	chomp;
	my @array=split "\t", $_;

	my $var_type=$array[0];
	my @var_list=split /,/, $array[1];
	my $list_len=$#var_list+1;

	print "Variable Type: '$array[0]'\n";
	print "List Length: $list_len\n";
	print "\n\n";

	if($var_type=~/groups/){
		my ($token, $gpname)=split /\//, $var_type;
		write_list_to_file(\@var_list, "$grp_dir/$gpname");
	}else{
		write_list_to_file(\@var_list, "$model_dir/$array[0]");
	}
	

}
close(MOD_FH);


###############################################################################

print STDERR "Done.\n";
