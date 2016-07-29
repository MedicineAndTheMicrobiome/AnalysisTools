#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_c $opt_f $opt_F $opt_o $opt_s $opt_p);
use Statistics::Descriptive;

my $SUBSET_SIZE=50;
my $FLOOR_QV = 20;

getopts("cf:F:o:s:p:");
my $usage = "usage: 
$0 
	-o <output file name>
	[-c input is CA frg file, instead of default fasta]
	[-f output a count of bases >= this floor value, default is 20]
	[-F output file name for bases >= the specified floor value]
	[-s <subset sample size, default $SUBSET_SIZE>]
        [-p <preference for summary of all but not the reads>]

	This program will read in a quality FASTA, or CA FRG file through STDIN and then
	output some some statistics on the sequences based on the quality.

	Output format is:

id\\tcount\\t[fullmean (stdev) min max]\\t[front_mean (stdev)]\\t[mid_mean (stdev)]\\t[end_mean (stdev)]\\n;

	Last line will be:

ALL\\t[Num_of_sequences]\\t[mean_length]\\t[fullmean (stdev min max)]\\t[front_mean (stdev)]\\t[mid_mean (stdev)]\\t[end_mean (stdev)]\\n;
";

if(!(defined($opt_o))){
	die $usage;
}

my $subset_size=$SUBSET_SIZE;
if(defined($opt_s)){
	$subset_size+=$opt_s;
}
if(defined $opt_f){
	$FLOOR_QV = $opt_f;
}
my $floor_fh;
if(defined $opt_F){
    open($floor_fh, ">$opt_F") || die "Could not open $opt_F for writing.\n";
}


# Global statistics
my $num_sequences=0;
my $all_lengths=0;
my $all_mean=0;
my $all_stdev=0;
my $allf_mean=0;
my $allf_stdev=0;
my $allm_mean=0;
my $allm_stdev=0;
my $alle_mean=0;
my $alle_stdev=0;
my ($allmin,$allmax) = (100000,0);
my $allfloor=0;

open(OUTPUT_FH, ">$opt_o") || die "Could not open $opt_o for writing.\n";

if (defined $opt_c) {
    processFrg();
} else {
    processFasta();
}

###############################################################################

sub processFrg {

    print STDERR "Processing FRG file...\n";

    my ($b,$e,$id);
    my @qvList;
    while(<>) {
        chomp;
        my $tag = substr($_,0,4);
        if($tag eq 'acc:')
        {
            $id=substr($_,4);
            @qvList = ();
        }
        elsif ( $tag eq 'clr:' ) {
        # clr comes last, after qlt, so process here
            ($b,$e) = split /,/,substr($_,4);
            my @clr = @qvList[$b..$e-1]; # space based coordinate conversion
            process_record( $id, \@clr );
        }
        elsif ( $tag eq 'qlt:' )
        {
            while(<>) {
                chomp;
                last if $_ eq '.';
                for my $qv (split '') {
                    push(@qvList, ord($qv) - ord(0));
                }
            }
        }
    }
}

###############################################################################

sub processFasta {

    print STDERR "Processing FASTA file...\n";

    my ($defline, $prev_defline, $sequence);
    while(<STDIN>){
        chomp;

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

###############################################################################
# Compute global statistics

# finish computation of averages
$all_lengths/=$num_sequences;
$all_mean/=$num_sequences;
$all_stdev/=$num_sequences;
$allf_mean/=$num_sequences;
$allf_stdev/=$num_sequences;
$allm_mean/=$num_sequences;
$allm_stdev/=$num_sequences;
$alle_mean/=$num_sequences;
$alle_stdev/=$num_sequences;

# trucante unnecessary digits
$all_lengths=sprintf("%5.2f",$all_lengths);
$all_mean=sprintf("%5.2f",$all_mean);
$all_stdev=sprintf("%5.2f",$all_stdev);
$allf_mean=sprintf("%5.2f",$allf_mean);
$allf_stdev=sprintf("%5.2f",$allf_stdev);
$allm_mean=sprintf("%5.2f",$allm_mean);
$allm_stdev=sprintf("%5.2f",$allm_stdev);
$alle_mean=sprintf("%5.2f",$alle_mean);
$alle_stdev=sprintf("%5.2f",$alle_stdev);

# output results
my $outstr="ALL\t$num_sequences\t$all_lengths\t[$all_mean ($all_stdev) $allmin $allmax]\t[$allf_mean ($allf_stdev)]\t[$allm_mean ($allm_stdev)]\t[$alle_mean ($alle_stdev)]\n";
print OUTPUT_FH $outstr;
close OUTPUT_FH;

if ( defined $floor_fh ) {
    print $floor_fh "ALL $allfloor\n";
    close $floor_fh;
}

print STDERR "Completed.\n";


###############################################################################
###############################################################################

sub getStats{
	my $values_ref=shift;
	my $stat=Statistics::Descriptive::Full->new();
	$stat->add_data(@{$values_ref});
	my $sum=$stat->sum();
	my $count=sprintf("%5i", $stat->count());
	my $mean=sprintf("%5.2f", $stat->mean());
	my $stdev=sprintf("%5.2f", $stat->standard_deviation());
	return($sum, $count, $mean, $stdev);
}

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

	# Compute region boundaries
	my $length=@$qv;
	my $front_end = $subset_size;
	my $middle_start = ($length -$subset_size)/2;
	my $middle_end = $middle_start+$subset_size;
	my $end_start = $length-$subset_size;
    my ($min,$max) = (100000,0);
    my $floor = 0;

	my (@front, @middle, @end);

	# Seperate regions into different datasets
	for(my $i=0; $i<@$qv; $i++){

        my $q = $qv->[$i];
        $max = $q if $q > $max;
        $min = $q if $q < $min;
        $floor++  if $q >= $FLOOR_QV;

		if($i<$front_end){
			push @front, $q;
		}	
		if($i>=$end_start){
			push @end, $q;
		}
		if($i>=$middle_start && $i<$middle_end){
			push @middle, $q;
		}
	}
    $allmax = $max if $max > $allmax;
    $allmin = $min if $min < $allmin;

    $allfloor += $floor;
    print $floor_fh "$id $floor\n" if defined $floor_fh;

	# Compute statistics for all regions
	my($sum, $count, $mean, $stdev)=getStats($qv);
	my($fsum, $fcount, $fmean, $fstdev)=getStats(\@front);
	my($msum, $mcount, $mmean, $mstdev)=getStats(\@middle);
	my($esum, $ecount, $emean, $estdev)=getStats(\@end);

	# Output per read statistics
    if (!defined $opt_p) {
	my $outstr="$id\t$length\t[$mean ($stdev) $min $max]\t[$fmean ($fstdev)]\t[$mmean ($mstdev)]\t[$emean ($estdev)]\n";
	print OUTPUT_FH $outstr;	
    } 
    


	# Sum up for global statistics
	$num_sequences++;
	$all_lengths+=$length;
	$all_mean+=$mean,
	$all_stdev+=$stdev;
	$allf_mean+=$fmean;
	$allf_stdev+=$fstdev;
	$allm_mean+=$mmean;
	$allm_stdev+=$mstdev;
	$alle_mean+=$emean;
	$alle_stdev+=$estdev;
}

#------------------------------------------------------------------------------
