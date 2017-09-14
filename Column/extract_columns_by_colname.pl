#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_k $opt_r $opt_K $opt_R $opt_d);

getopts(":i:k:r:K:R:d:");

my $DELIM="\t";

my $usage = "
	usage:
	$0
		-i <input table file>

		One of:
		(File)
		[-k <keep target column names, filename>]
		[-r <remove target column names, filename>]
		(List)
		[-K <list of colnames to Keep, comma separated (command-line)>]
		[-R <list of colnames to Remove, comma separated (command-line)>]

		[-d <column delimiter, default=tab>]
		
	Output goes into STDOUT.

	You can put comments in our keep/remove file.  Anything after and including a # will be 
	considered comment.

	For example:

	Variable1
	# Ignore
	Variable2 # Ignore

	Is the same was:
	Variable1
	Variable2
	
";

if(!defined($opt_i)){
	die $usage;
}

my $InFile=$opt_i;

my $keep=0;
my $file="";
my $colnames_ref;

###############################################################################

sub read_file{
	my $file=shift;
	my @list;
	open(FH, "<$file") || die "Could not open $file for reading.\n";
	while(<FH>){
		chomp;

		# Remove comments
		my @values=split "#", $_;
		
		# Remove leading/trailing spaces
		my $item=$values[0];
		$item=~s/^\s+//;
		$item=~s/\s+$//;

		push @list, $item;
	}
	close(FH);
	return(\@list);
}

sub get_list{
	my $string=shift;
	my @arr=split /,\s*/, $string;
	return(\@arr);
}

###############################################################################

if(defined($opt_k)){
	$keep=1;
	$colnames_ref=read_file($opt_k);
}
if(defined($opt_r)){
	$keep=0;
	$colnames_ref=read_file($opt_r);
}

if(defined($opt_K)){
	$keep=1;
	$colnames_ref=get_list($opt_K);
}
if(defined($opt_R)){
	$keep=0;
	$colnames_ref=get_list($opt_R);
}

if(defined($opt_d)){
	$DELIM=$opt_d;
}
###############################################################################

sub print_arr{
	my $arr_ref=shift;
	foreach my $val (@{$arr_ref}){
		print STDERR "--> $val\n";
	}
}

###############################################################################

if($keep){
	print STDERR "Keeping:\n";
}else{
	print STDERR "Removing:\n";
}

print_arr($colnames_ref);
my %colname_hash;
foreach my $cn(@{$colnames_ref}){
	$colname_hash{$cn}=1;	
}

###############################################################################

open(IN_FH, "<$InFile") || die "Could not open $InFile\n";

$_=<IN_FH>;
chomp $_;
my @header_arr=split /$DELIM/, $_;
my $num_col=$#header_arr+1;

# Determine which columsn to keep:
my @keep_ix;
for(my $i=0; $i<$num_col; $i++){
	my $hdr=$header_arr[$i];
	if($colname_hash{$hdr}){
		if($keep==1){
			push @keep_ix, $i;			
		}
	}else{
		if($keep==0){
			push @keep_ix, $i;
		}
	}
}

my $num_kept=$#keep_ix +1;
print STDERR "Number of Columns to Keep: $num_kept\n";

# Output header
my @out_arr;
foreach my $i (@keep_ix){
	push @out_arr, $header_arr[$i];
}
print STDOUT (join $DELIM, @out_arr) . "\n";

###############################################################################

# Output values
while(<IN_FH>){

	chomp;
	my @cols_arr=split /$DELIM/, $_;

	my @out_arr;
	foreach my $i (@keep_ix){
		push @out_arr, $cols_arr[$i];
	}
	
	print STDOUT (join $DELIM, @out_arr) . "\n";

}

close(IN_FH);

###############################################################################

print STDERR "Done.\n";
