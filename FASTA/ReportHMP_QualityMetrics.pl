#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_r $opt_c $opt_b $opt_o $opt_p $opt_q);
use Statistics::Descriptive;

my $MIN_READ_LEN=300;
my $MIN_20_QV = 300;
my $MIN_READS = 3000;

getopts("b:p:o:r:c:q:");
my $usage = "usage:
$0
        -p <project>
        -o <output file name>
        -b <barcode sample association file>
        [-r input is minimum reads in a file, instead of default 3000]
        [-c input a count of bases in a read>, default is 300]
        [-q input is minimum 20 qv bases in a read, default is 300]


        This program will read in a quality FASTA file through STDIN and then
        output some some statistics on the sequences based on the quality.

        Output format is:

id\\t total count \\t good reads \\t reads less than min length \\t reads with min Q20 bases\\n

        Last line will be:

Number of Samples less than MIN_NUM_READS  are MIN_READ_COUNT out of  TOTAL samples\n\n
";

if(!(
	defined($opt_o) &&
	defined($opt_b) &&
	defined($opt_p)
   )){
        die $usage;
}

if(defined($opt_c)){
    $MIN_READ_LEN = $opt_c;
}
if(defined $opt_q){
    $MIN_20_QV = $opt_q;
}
if(defined $opt_r){
    $MIN_READS = $opt_r;
}
my $name=$opt_p;


# Load barcodes
print STDERR "Loading Barcodes...\n";
my $barcode_hash_ref=load_barcodes($opt_b);
my %barcode_hash=%$barcode_hash_ref;
my @barcodes_arr=sort keys %{$barcode_hash_ref};
my $num_barcodes=$#barcodes_arr+1;
my $min_read_length = 0;
my $discrad_seq = 0;
my $num_sequences = 0;
my $min_read_count = 0;

my $outListFile=$opt_o;
$outListFile =~s/qual\.summary//;
my $outListFile_name=$outListFile."Failed.list";
my (@Failed_array, @noamplicon_array, $noampliconCount, $flag);

print STDERR "Num barcodes loaded: $num_barcodes\n";

open(OUTPUT_FH, ">$opt_o") || die "Could not open $opt_o for writing.\n";

foreach my $bar (keys %barcode_hash){
    my $input_file = "$name\.$barcode_hash{$bar}\.fasta\.qual";
    $min_read_length = 0;
    $discrad_seq = 0;
    $num_sequences = 0;
    my $num_bad_seq = 0;
    my $flag = 0;
    my $num_good_seq = 0;
    if (-e $input_file){
	my $defline;
	my $sequence;
	my $prev_defline;
	
	
	open(IN_FH, "<$input_file") || die "Could not open $input_file\n";
	while(<IN_FH>){
	    chomp;
	    # Initialize Global Statistics
	   
	    if(/^>/){
		$defline=$_;
		if($sequence ne ""){
		    processFastaRecord($prev_defline, $sequence);
		    $sequence="";
		}
		$prev_defline=$defline;
	    }else{
		$sequence.=($_ . " ");
	    }
	}
	processFastaRecord($prev_defline, $sequence);
    }
    else{
	print STDERR "The file $input_file does not exist\n";
	push @noamplicon_array,$barcode_hash{$bar};
    }
# Compute Global Satistics
    
    $num_bad_seq = $min_read_length + $discrad_seq;
    $num_good_seq = $num_sequences - $num_bad_seq;

    my $out_result = "$barcode_hash{$bar}\t$num_sequences\t$num_good_seq\t$min_read_length\t$discrad_seq\n";
    

    print OUTPUT_FH $out_result;
    if($num_good_seq < $MIN_READS){
	$min_read_count++;
	$flag=1;
    }

    if ($flag == 1){
	push @Failed_array,$barcode_hash{$bar}; 
    }

    close(IN_FH);
}

$noampliconCount=$#noamplicon_array+1;
print OUTPUT_FH "Number of Samples less than $MIN_READS reads are $min_read_count and $noampliconCount had no amplicon out of $num_barcodes samples\n";
close(OUTPUT_FH);

open (OUT,">$outListFile_name")|| die "Could not open $outListFile_name for writing.\n";

print OUT "Samples which had no amplicons in run\n";

foreach my $list(@noamplicon_array){
    print OUT "$list\n";
}

print OUT "Samples found but failed the HMP Metrics\n";

foreach my $list(@Failed_array){
    print OUT "$list\n";
}

close(OUT);

###################################################################################################
sub load_barcodes{
        # Loads barcodes and sample from file into hash
        my $filename=shift;
        my %hash;
        open(FH, "<$filename") || die "Could not open $filename\n";
        while(<FH>){
                chomp;
                my ($barcode, $id)=split /\t/, $_;
                $hash{$barcode}=$id;
                }
        return(\%hash);
}
# End sub load_barcodes

###################################################################################################
sub processFastaRecord($$) {
        my $defline = shift;
        my $sequence = shift;

        # Parse out id from defline
        my $id;
        if($defline=~/^>(\S+)/){
                $id=$1;
        }

        # Convert quality value string into array
        $sequence=~s/^\s+//;
        my @qv=split /\s+/, $sequence;
    process_record($id,\@qv);
}
sub process_record($\@){
    my $id = shift;
    my $qv = shift;
    
    # Compute the quals
    my $length = @$qv;
    if ($length < $MIN_READ_LEN){
	$min_read_length++;
    }
    else{
	my $count_20_qual = 0;
	for(my $i=0; $i<@$qv; $i++){
	    my $q = $qv->[$i];
	    if ($q < "20"){
		$count_20_qual++;
	    }
	}
	if($count_20_qual > $MIN_20_QV){
	    $discrad_seq++;
	}
    }
	
    $num_sequences++;
}
##########################################################################################################
