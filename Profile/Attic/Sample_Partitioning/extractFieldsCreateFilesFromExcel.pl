#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_l);
getopts("i:l:");

my $usage = "usage:
$0
        Find fields and generate files for those fields.
                -i <input tab delimited sample attribute file>
                -l <comma separated list of columns eg. region,type>
        NOTE: The first column in input file should be sample Id
        Gives a list of samples in each file combination.

";

if(!(
        defined($opt_i) &&
	defined($opt_l)
)){
        die $usage;
}

my $num=0;
my $flag=0;
my $input_file_name=$opt_i;
my @strings=split (',', $opt_l);
my %col_found;
my (%ld_found,%field_found,$id,@final_array, %lds);
print STDERR "Input Filename: $input_file_name\n";


open(IN_FH, "<$input_file_name") || die "Could not open $input_file_name\n";

while(<IN_FH>){
    chomp;
    my @fields = split /\t/, $_;
    my $j=0;
    my @arr;
    my $combination;
    if ($num==0){
	foreach my $str(@strings){
	    for(my $i=0; $i<=$#fields; $i++){
		if ($fields[$i]=~/$str/i){
		    $col_found{$str}=$i;
		}
		
	    }
	}
    }

    else{
	foreach my $key(keys %col_found){
	    my @id;
	    $id=$fields[$col_found{$key}];
	    if($id =~ /\//){
		@id=split('/', $id);
		$id=join '-',@id;
	    }
	    else{
		@id=split(' ', $id);
		$id=join '-',@id;
	    }
	    push (@arr,$id);
	    $field_found{$key}=$id;
	}
	$combination=join '_',@arr;
	push (@final_array,$combination);
	$ld_found{$combination}=1;
	} 
    if(defined(	$ld_found{$combination})){
	push (@{$lds{$combination}},$fields[0]); 
       }
 $num++;
}
print "Total number of rows scanned = $num\n";
close(IN_FH);

my $basename=$input_file_name;
my @out_names;
my @filehandles;
my $i=0;
my @id_array=keys(%ld_found);
my $i=0;

foreach my $split_attribute(sort @id_array){
    my $out_names=$basename . ".".$split_attribute.".list";
    print "Output Files: '$out_names'\n";
    open(OUT_FH,">$out_names") || die "can't open '$out_names': [$!]\n";
    foreach(@{$lds{$split_attribute}}){
	print OUT_FH $_,"\n";
    }
    close(OUT_FH);
}



		
	    


