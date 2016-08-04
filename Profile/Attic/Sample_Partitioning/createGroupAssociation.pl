#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_g $opt_a $opt_s);
getopts("i:g:a:s:");

my $usage = "usage:
$0
       Takes the following input
                -i <input tab delimited sample attribute file>
                -g <group name>
                -a <list of comma separated associations with group>
                -s <comma separated special column name,attribute of column>
        NOTE: The first column in input file should be sample Id
       
        Gives a comma separated list of samples for each group association.
        Writes out:
                <input file><group>_association.list

";

if(!(
        defined($opt_i) &&
	defined($opt_g) &&
	defined($opt_a)
)){
        die $usage;
}

my $input_file_name=$opt_i;

print "Input file :$input_file_name\n";

my $group= $opt_g;
my @strings=split (',', $opt_a);
my (@special,$spcl_attri, $spcl_col, $spcl_col_found);

if(defined($opt_s)){
  @special=split(',',$opt_s);
  $spcl_col=$special[0];
  $spcl_attri=$special[1];
}

my $num=0;
my $j=0;
my %grp_hash;
my ($combination, $grp_id, $id_spcl);
my (%ld_found,%field_found,$id,@final_array, %lds,%grp_asso_hash,%spcl_combination_hash, %spcl_ld_found);
my (%col_found,$grp_found,@assoc_array,%combination_hash, $spcl_grp_id);
open(IN_FH, "<$input_file_name") || die "Could not open $input_file_name\n";

while(<IN_FH>){
    chomp;
    my @arr;
    my @fields = split /\t/, $_;
    if ($num==0){
	foreach my $str(@strings){
	    for(my $i=0; $i<=$#fields; $i++){
		$col_found{$str}=$i if ($fields[$i]=~/$str/i);
		$grp_found=$i if($fields[$i]=~/$group/i);
		$spcl_col_found=$i if($fields[$i]=~/$spcl_col/i);
	    }
	}
    }
    else{
	$id_spcl=$fields[$spcl_col_found];
	$grp_id=$fields[$grp_found];
	if ($grp_id ne ''){
	    $grp_hash{$grp_id}=1;
	    foreach my $key(keys %col_found){
		$id=$fields[$col_found{$key}];
		my @id=split(' ', $id);
		$id=join '-',@id;
		push (@arr,$id);
	    }
	    $combination=join '_',@arr;
	    $combination_hash{$combination}=1;
	    if ((defined($opt_s)) && ($id_spcl eq $spcl_attri)){
		#$spcl_grp_id = $grp_id . "_" . $id_spcl;
		$spcl_ld_found{$grp_id}{$combination}=$fields[0];
	    }
	    else {
		$ld_found{$grp_id}{$combination}=$fields[0];
	    }		
	}
    }
    $num++;
}
close(IN_FH);


my $out_file=$input_file_name.".".$group."_association.list";

print "Output File Name : $out_file\n";

open(OUT_FH,">$out_file") || die "can't open '$out_file': [$!]\n";

print OUT_FH "GroupId";

if(defined($opt_s)){
    foreach my $newkey(sort keys %combination_hash){   
	print OUT_FH ",$newkey";
    }
    print OUT_FH "\n";

    foreach my $key (sort keys %spcl_ld_found){
	print OUT_FH "$key";
	foreach my $newkey(sort keys %combination_hash){
	    if (defined($spcl_ld_found{$key}{$newkey})){
		print OUT_FH ",$spcl_ld_found{$key}{$newkey}";
	    }
	    else{
		print OUT_FH ",";
	    } 
	}
	    print OUT_FH "\n";
    }   
}
else{
    foreach my $newkey(sort keys %combination_hash){
	print OUT_FH ",$newkey";
    }

    print OUT_FH "\n";
    
    foreach my $key (sort keys %ld_found){
	print OUT_FH "$key";
	foreach my $newkey(sort keys %combination_hash){
	    if (defined($ld_found{$key}{$newkey})){
		print OUT_FH ",$ld_found{$key}{$newkey}";
	    }
	    else{
		print OUT_FH ",";
	    } 
	}
	print OUT_FH "\n";
    }
}

	    
