#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use Statistics::Distributions;
use vars qw($opt_a $opt_b $opt_o);
getopts("a:b:o:");

my @filter_levels=(0, 0.0125, 0.025, .05, .1, .25, .333, .5);
my $filter_string=join ", ", @filter_levels;

my $usage = "usage:
$0
	-a <sample set A summary table>
	-b <sample set B summary table>
	-o <output>

	Reads in two summary tables, and determines what is in:

	Only A, Only B, A&B, and (A|B)'

	This will be computed for counts after filtering at percent levels:
	
	$filter_string

	The output will be two files, <output>.venn and <output>.stat
	Both outputs can be opened with Excel, or openoffice.org.

	The .venn file contains the categories between A and B that are 
	shared or exclusive depending on the cutoff.  The columns are:
	<category>\\t<A count>\\t<P[A]>\\t<B count>\\t<P[B]>\\t<Z-test statistic>\\t<P-value>\\n

	The .stat file contains the statistical significance of the differences between
	A and B, sorted by the Z-test statistic.
	<Category>\\t<A Count>\\t<P[A]>\\t<B Count>\\t<P[B]>\\t<A/B Ratio>\\t<A-B Diff>\\t<Z-test>\\t<Stat>\\t<P-value>\\t<Flag>\\n

	If the test statistic is less than 0, then P[A]<P[B].
	If the test statistic is greater than 0, then P[A]>P[B].
	The Flag notes the following conditions:
		*** = p-value <= .001
		**  = p-value <= .01
		*   = p-value <= .05

	The smaller the p-value (more extreme the Z-test statistic),
	the more likely the two distributions for the category are not the same.


";

if(!(
	defined($opt_a) &&
	defined($opt_b) &&
	defined($opt_o)
)){
	die $usage;
}


my $fileA=$opt_a;
my $fileB=$opt_b;
my $outfile=$opt_o;

my ($Acat_to_count_hash_ref)=load_samples($fileA);
my ($Bcat_to_count_hash_ref)=load_samples($fileB);

my ($normalAref, $Atotal)=normalize($Acat_to_count_hash_ref);
my ($normalBref, $Btotal)=normalize($Bcat_to_count_hash_ref);

my ($cat_to_test_statistic_hash_ref, $cat_to_pvalue_hash_ref)=
	compute_test_stat($Atotal, $Btotal, $normalAref, $normalBref, $Acat_to_count_hash_ref, $Bcat_to_count_hash_ref);

my $numAcat=scalar(keys %{$Acat_to_count_hash_ref});
my $numBcat=scalar(keys %{$Bcat_to_count_hash_ref});

#for(my $i=0; $i<=$#{$catAref}; $i++){
#	my $norm_count_to_show=sprintf("%3.5f", ${$normalAref}[$i]);
#	print "${$catAref}[$i] \t ${$countAref}[$i] $norm_count_to_show\n";
#}

#--------------------------------------------------------------------------

my $fh = FileHandle->new;
$fh->open(">$outfile\.venn") || die "Could not open $outfile\.venn\n";

print {$fh} "Num A Categories:\t $numAcat\n";
print {$fh} "Num B Categories:\t $numBcat\n";
print {$fh} "\n\n";

foreach my $level (@filter_levels){

	my ($Aonly, $Bonly, $Both, $Neither)=compare($normalAref, $normalBref, $level);

	output_count_venn_results(
		$fh,
		$level,
		$Aonly, $Bonly, $Both, $Neither, 
		$Atotal, $Btotal, 
		$Acat_to_count_hash_ref, $Bcat_to_count_hash_ref,
		$cat_to_test_statistic_hash_ref, $cat_to_pvalue_hash_ref,
		$fileA, $fileB
		);
	
	
	my $aonlycount=$#{$Aonly}+1;
	my $bonlycount=$#{$Bonly}+1;
	my $bothcount=$#{$Both}+1;
	my ($jaccard, $sorenson, $kulczynski)=compute_similarity_indices($aonlycount+$bothcount, $bonlycount+$bothcount, $bothcount);

	print {$fh} "\n";
	print {$fh} "Jaccard coef =\t$jaccard\n";
	print {$fh} "Sorensen's coef =\t$sorenson\n";
	print {$fh} "Kulczynski's coef =\t$kulczynski\n";
	print {$fh} "\n\n\n";
	
}

$fh->close;

#--------------------------------------------------------------------------

my $fh = FileHandle->new;
$fh->open(">$outfile\.stat") || die "Could not open $outfile\.stat\n";

output_statistically_significant_differences(
	$fh,
	$Acat_to_count_hash_ref, $Bcat_to_count_hash_ref,
	$normalAref, $normalBref,
	$cat_to_test_statistic_hash_ref,
	$cat_to_pvalue_hash_ref,
	$fileA, $fileB	
);

$fh->close;

print STDERR "Done.\n";

###########################################################################

sub compute_similarity_indices{
	my $num_a=shift;
	my $num_b=shift;
	my $num_aAndb=shift;
	
	my $jaccard=compute_jaccard($num_a, $num_b, $num_aAndb);
	my $sorenson=compute_sorenson($num_a, $num_b, $num_aAndb);
	my $kulczynski=compute_kulczynski($num_a, $num_b, $num_aAndb);

	return($jaccard, $sorenson, $kulczynski);

}

sub compute_jaccard{
	my $a=shift;
	my $b=shift;
	my $ab=shift;
	if(($a+$b-$ab)==0){
		return("Undef");
	}else{
		return($ab/($a+$b-$ab));
	}
}

sub compute_sorenson{
	my $a=shift;
	my $b=shift;
	my $ab=shift;
	if(($a+$b)==0){
		return("Undef");
	}else{
		return(2*$ab/($a+$b));
	}
}

sub compute_kulczynski{
	my $a=shift;
	my $b=shift;
	my $ab=shift;
	if($a==0 || $b==0){
		return("Undef");
	}else{
		return((($ab/$a) + ($ab/$b))/2.0);
	}
}


sub print_hash{
	my $hash_ref=shift;
	foreach my $id(sort keys %{$hash_ref}){
		print "$id: ${$hash_ref}{$id}\n";
	}
}

sub output_statistically_significant_differences{
	my $fh=shift;
	my $a_cat2count=shift;
	my $b_cat2count=shift;
	my $a_norm_hash_ref=shift;
	my $b_norm_hash_ref=shift;
	my $test_stat_hash_ref=shift;
	my $pvalue_hash_ref=shift;
	my $fileA=shift;
	my $fileB=shift;

	my %anorm_hash;
	my %bnorm_hash;

	my @info_arr;
	foreach my $cat (keys %{$test_stat_hash_ref}){

		my $acat_count=${$a_cat2count}{$cat};
		my $bcat_count=${$b_cat2count}{$cat};

		my $Aperc=${$a_norm_hash_ref}{$cat};
		my $Bperc=${$b_norm_hash_ref}{$cat};

		# Compute ratio of A to B
		my $ratio;
		if($bcat_count>0){
			$ratio=$Aperc/$Bperc;
		}else{
			$ratio="Inf";
		}

		# Compute flag based on p-value
		my $pvalue=${$pvalue_hash_ref}{$cat};
		my $flag="";
		if(defined($pvalue)){
			if($pvalue <= .001){
				$flag="***";
			}elsif($pvalue <= .01){
				$flag="**";
			}elsif($pvalue <= .05){
				$flag="*";
			}
		}

		# Store column info into a single line
		my @columns=($cat, $acat_count, $Aperc, $bcat_count, $Bperc, $ratio, ($Aperc-$Bperc), ${$test_stat_hash_ref}{$cat}, ${$pvalue_hash_ref}{$cat}, $flag);
		push @info_arr, \@columns;
	}

	# Sort by test statistic
	my @sorted_info_arr=sort by_test_stat @info_arr;

	# Output results
	print {$fh} "Category\tA Count\tP[A]\tB Count\tP[B]\tA/B Ratio\tA-B Diff\tZ-test Stat\tP-value\tFlag\n";
	foreach my $info(@sorted_info_arr){
		print {$fh} (join "\t", @{$info}) . "\n";
	}

	print {$fh} "\n\n";
	print {$fh} "*** = p-value <= .001\n";
	print {$fh} "**  = p-value <= .01\n";
	print {$fh} "*   = p-value <= .05\n";

}

sub by_test_stat{
	${$a}[6] <=> ${$b}[6];
}

###########################################################################

sub output_count_venn_results{
	my $fh=shift;
	my $cutoff=shift;
	my $aonly=shift;
	my $bonly=shift;
	my $both=shift;
	my $neither=shift;
	my $atotal=shift;
	my $btotal=shift;
	my $a_cat2count=shift;
	my $b_cat2count=shift;
	my $test_stat_hash_ref=shift;
	my $pvalue_hash_ref=shift;
	my $fileA=shift;
	my $fileB=shift;

	my $aonly_count=$#{$aonly}+1;
	my $bonly_count=$#{$bonly}+1;
	my $both_count=$#{$both}+1;
	my $neither_count=$#{$neither}+1;

	my %anorm_hash;
	my %bnorm_hash;

	foreach my $cat(keys %{$a_cat2count}){
		$anorm_hash{$cat}=sprintf("%1.5f", ${$a_cat2count}{$cat}/$atotal);
	}
	foreach my $cat(keys %{$b_cat2count}){
		$bnorm_hash{$cat}=sprintf("%1.5f", ${$b_cat2count}{$cat}/$btotal);
	}

	print {$fh} "Cutoff:\t$cutoff\tA Count\tP[A]\tB Count\tP[B]\tZ-Test Stat\tP-Value\n";
	print {$fh} "$fileA [$aonly_count]\n";
	for(my $i=0; $i<=$#{$aonly}; $i++){
		my $cat=${$aonly}[$i];
		my $acat_count=${$a_cat2count}{$cat};
		my $bcat_count=${$b_cat2count}{$cat};
		my $Aperc=$anorm_hash{$cat};
		my $Bperc=$bnorm_hash{$cat};
		print {$fh} "\t$cat\t$acat_count\t$Aperc\t$bcat_count\t$Bperc\t${$test_stat_hash_ref}{$cat}\t${$pvalue_hash_ref}{$cat}\n";
	}

	print {$fh} "$fileB [$bonly_count]\n";
	for(my $i=0; $i<=$#{$bonly}; $i++){
		my $cat=${$bonly}[$i];
		my $acat_count=${$a_cat2count}{$cat};
		my $bcat_count=${$b_cat2count}{$cat};
		my $Aperc=$anorm_hash{$cat};
		my $Bperc=$bnorm_hash{$cat};
		print {$fh} "\t$cat\t$acat_count\t$Aperc\t$bcat_count\t$Bperc\t${$test_stat_hash_ref}{$cat}\t${$pvalue_hash_ref}{$cat}\n";
	}

	print {$fh} "$fileA/$fileB [$both_count]\n";
	for(my $i=0; $i<=$#{$both}; $i++){
		my $cat=${$both}[$i];
		my $acat_count=${$a_cat2count}{$cat};
		my $bcat_count=${$b_cat2count}{$cat};
		my $Aperc=$anorm_hash{$cat};
		my $Bperc=$bnorm_hash{$cat};
		print {$fh} "\t$cat\t$acat_count\t$Aperc\t$bcat_count\t$Bperc\t${$test_stat_hash_ref}{$cat}\t${$pvalue_hash_ref}{$cat}\n";
	}

	print {$fh} "Neither [$neither_count]\n";
	for(my $i=0; $i<=$#{$neither}; $i++){
		my $cat=${$neither}[$i];
		my $acat_count=${$a_cat2count}{$cat};
		my $bcat_count=${$b_cat2count}{$cat};
		my $Aperc=$anorm_hash{$cat};
		my $Bperc=$bnorm_hash{$cat};
		print {$fh} "\t$cat\t$acat_count\t$Aperc\t$bcat_count\t$Bperc\t${$test_stat_hash_ref}{$cat}\t${$pvalue_hash_ref}{$cat}\n";
	}
	
}

###########################################################################

sub compute_test_stat{
	my $Atotal=shift;
	my $Btotal=shift;
	my $normAref=shift;
	my $normBref=shift;
	my $cat_to_countAref=shift;
	my $cat_to_countBref=shift;

	my %a_normal=%{$normAref};
	my %b_normal=%{$normBref};

	# Put all categories in the same hash so we know we can cycle through all the categories
	my %all_cat;
	foreach my $cat(keys %{$cat_to_countAref}){
		$all_cat{$cat}=1;
	}
	foreach my $cat(keys %{$cat_to_countBref}){
		$all_cat{$cat}=1;
	}

	my %cat_to_teststat_hash;
	my %cat_to_pvalue_hash;
	foreach my $cat(keys %all_cat){
		my $phat=(${$cat_to_countAref}{$cat}+${$cat_to_countBref}{$cat})/($Atotal+$Btotal);
		my $num;
		my $den;
		if($phat!=0){
			$cat_to_teststat_hash{$cat}=($a_normal{$cat}-$b_normal{$cat})/sqrt($phat*(1-$phat)*((1/$Atotal) + (1/$Btotal)));
			$cat_to_pvalue_hash{$cat}=2*(Statistics::Distributions::uprob(abs($cat_to_teststat_hash{$cat})));
		}else{
			$cat_to_teststat_hash{$cat}=undef;
			$cat_to_pvalue_hash{$cat}=undef;
		}
	}

	return(\%cat_to_teststat_hash, \%cat_to_pvalue_hash);
}

###########################################################################

sub compare{
	my $normAref=shift;
	my $normBref=shift;
	my $cutoff=shift;
	
	my @Aonly;
	my @Bonly;
	my @Both;
	my @Neither;

	my %a_candidates;
	my %b_candidates;
	my %all;

	foreach my $cat(keys %{$normAref}){
		if(${$normAref}{$cat}>$cutoff){
			$a_candidates{$cat}++;	
		}
		$all{$cat}=1
	}

	foreach my $cat(keys %{$normBref}){
		if(${$normBref}{$cat}>$cutoff){
			$b_candidates{$cat}++;	
		}
		$all{$cat}=1
	}

	foreach my $candidate(keys %all){
		if($a_candidates{$candidate} && $b_candidates{$candidate}){
			push @Both, $candidate;
		}elsif($a_candidates{$candidate}){
			push @Aonly, $candidate;
		}elsif($b_candidates{$candidate}){
			push @Bonly, $candidate;	
		}else{
			push @Neither, $candidate;
		}
	}
	
	@Aonly=sort @Aonly;
	@Bonly=sort @Bonly;
	@Both=sort @Both;
	@Neither=sort @Neither;
	
	return(\@Aonly, \@Bonly, \@Both, \@Neither);
}

###########################################################################

sub normalize{
	my $hash_ref=shift;
	my %normalized_hash;
	
	my $sum=0;
	foreach my $cat(keys %{$hash_ref}){
		$sum+=${$hash_ref}{$cat};
	}

	foreach my $cat(keys %{$hash_ref}){
		$normalized_hash{$cat}=${$hash_ref}{$cat}/$sum;
	}
	return(\%normalized_hash, $sum);
}

###########################################################################

sub load_samples{
	my $filename=shift;
	open(FH, "<$filename") || die "Could not open $filename\n";

	my $first_line=<FH>;
	chomp $first_line;

	my @categories=split /\t/, $first_line;
	shift @categories;	# Get rid of "Sample"
	shift @categories;	# Get rid of "Total"

	my @total_counts;
	while(<FH>){
		chomp;
		my @counts=split /\t/, $_;
		for(my $i=2; $i<=$#counts; $i++){
			$total_counts[$i-2]+=$counts[$i]
		}
	}

	my %cat_to_count_hash;
	for(my $i=0; $i<=$#total_counts; $i++){
		$cat_to_count_hash{$categories[$i]}=$total_counts[$i];
	}

	return(\%cat_to_count_hash);
}


