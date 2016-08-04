#!/usr/bin/env perl

##########################################################################################

use strict;
use Getopt::Std;
use Statistics::Distributions;


use vars qw($opt_s $opt_c $opt_t $opt_o);
getopts("s:c:t:o:");

my $usage = "usage:
$0

	-s <summary_table.xls file, generated from Pull_Taxonomic_Crosssection_By_Rank.pl>
	-c <clusters file to compare within>
	-t <percentage threshold to compute p-values, for example 5, for 5%>  
	-o <output file>

";

if(!(
	defined($opt_s) &&
	defined($opt_c) &&
	defined($opt_t) &&
	defined($opt_o)
)){
	die $usage;
}

my $summary_table_fn=$opt_s;
my $compare_clust_fn=$opt_c;
my $percent=$opt_t;
my $output_fn=$opt_o;


my ($column_to_category_map_ref, $category_to_column_map_ref, $sample_counts_hash_ref)=load_summary_table($summary_table_fn);
my ($clusters_hash_ref)=load_clusters($compare_clust_fn);

foreach my $cluster(keys %{$clusters_hash_ref}){
	print STDERR "\n\nWorking on cluster: $cluster\n";
	run_tests(${$clusters_hash_ref}{$cluster}, $sample_counts_hash_ref);
}

##########################################################################################

sub run_tests{
	my $cluster_arr_ref=shift;
	my $sample_counts_arr_ref=shift;

	foreach my $member(@{$cluster_arr_ref}){
		print STDERR "\t$member\n";
	}

	print STDERR "\n\n Computing homogeneity across cluster\n";
	chisqr_test_of_homogeneity($cluster_arr_ref, $sample_counts_hash_ref);

	print STDERR "\n\n Computing Ratios\n";
	my ($versus_all_results, $versus_another_results)=chisqr_test_of_ratios($cluster_arr_ref, $sample_counts_hash_ref, $percent);

	foreach my $results(@{$versus_all_results}){
		print "$results\n";
	}

	foreach my $results(@{$versus_another_results}){
		print "$results\n";
	}

}

sub chisqr_test_of_ratios{
	my $cluster_members_arr_ref=shift;
	my $sample_counts_hash_ref=shift;
	my $threshold_cutoff=shift;
	$threshold_cutoff/=100.0;

	my @members=@{$cluster_members_arr_ref};	

	my $num_members=$#members+1;

	my @members_counts;
	foreach my $member(@members){
		print STDERR "input counts [$member]: " . join (",", @{${$sample_counts_hash_ref}{$member}}) . "\n";
		push @members_counts, ${$sample_counts_hash_ref}{$member};
	}
	my $num_categories=$#{$members_counts[0]}+1;
	
	print STDERR "num members: $num_members\n";
	print STDERR "num categories: $num_categories\n";

	my @totals;
	my $all_total=0;
	my @member_totals;
	my $member_id=0;
	foreach my $count_ref(@members_counts){
		for(my $cat=0; $cat<$num_categories; $cat++){
			$totals[$cat]+=${$count_ref}[$cat];
			$all_total+=${$count_ref}[$cat];
			$member_totals[$member_id]+=${$count_ref}[$cat];
		}	
		$member_id++;
	}

	print STDERR "total: \n\t" .  (join ",", @totals) . "\n";
	print STDERR "all total: $all_total\n";
	
	my @percentages;
	my @perc_forshow;
	for(my $cat=0; $cat<$num_categories; $cat++){
		$percentages[$cat]=$totals[$cat]/$all_total;
		$perc_forshow[$cat]=sprintf("%1.4f",$percentages[$cat]);
	}
	print STDERR "percs: \n\t" .  (join ",", @perc_forshow) . "\n";

	my @to_compute;
	for(my $cat=0; $cat<$num_categories; $cat++){
		$to_compute[$cat]=($percentages[$cat]>$threshold_cutoff)?1:0;
	}
	print STDERR "to compute: \n\t" .  (join ",", @to_compute) . "\n";

	print STDERR "Member totals: \n\t" .  (join ",", @member_totals) . "\n";

	
	# Compute ratios of specific categories to all
	print STDERR "Computing category X versus All\n";
	my @versus_all_results;
	for(my $cat=0; $cat<$num_categories; $cat++){
		if($to_compute[$cat]==1){

			my @member_ratios_arr;

			# Get observed ratios
			for(my $member_id=0; $member_id<$num_members; $member_id++){
				# Get ratio of category to all
				my $ratio_ref=get_ratio($members_counts[$member_id], $cat);
				print STDERR "obs ratio category [$cat]:\n\t" . (join ",",@{$ratio_ref}) . "\n";

				push @member_ratios_arr, $ratio_ref;
			}

			my ($chisqr_stat, $df, $pvalue)=compute_homogeneity_pvalue(\@member_ratios_arr);

			push @versus_all_results, (join ",",($cat, $chisqr_stat, $df, $pvalue));
			print STDERR "\n";
		}
	}

	# Compute ratios of specific categories pair wise
	print STDERR "Computing category X versus Y\n";
	my @versus_other_results;
	for(my $cat1=0; $cat1<$num_categories; $cat1++){
		for(my $cat2=0; $cat2<$num_categories; $cat2++){
			if(($to_compute[$cat1]==1) && ($to_compute[$cat2]==1) && ($cat1<$cat2)){
				
				my @member_ratios_arr;
				foreach(my $member_id=0; $member_id<$num_members; $member_id++){
					my @ratio_arr;
					push @ratio_arr, ${$members_counts[$member_id]}[$cat1];
					push @ratio_arr, ${$members_counts[$member_id]}[$cat2];
					push @member_ratios_arr, \@ratio_arr;
				}

				my ($chisqr_stat, $df, $pvalue)=compute_homogeneity_pvalue(\@member_ratios_arr);

				push @versus_other_results, (join ",",("$cat1/$cat2", $chisqr_stat, $df, $pvalue));
				print STDERR "\n";
			}
		}
	}

	return(\@versus_all_results, \@versus_other_results);
	

}

sub get_ratio{
	my $arr_ref=shift;
	my $target_idx=shift;
	
	my $sum=0;
	my $target;

	$target=${$arr_ref}[$target_idx];
	for(my $cat=0; $cat<=$#{$arr_ref}; $cat++){
		$sum+=${$arr_ref}[$cat];
	}

	my @ratio=($target, $sum-$target);
	return(\@ratio);
}


sub chisqr_test_of_homogeneity{
	my $cluster_members_arr_ref=shift;
	my $sample_counts_hash_ref=shift;

	my @members=@{$cluster_members_arr_ref};
	my @members_counts;
	foreach my $member(@members){
		print STDERR "input counts [$member]: " . join (",", @{${$sample_counts_hash_ref}{$member}}) . "\n";
		push @members_counts, ${$sample_counts_hash_ref}{$member};
	}

	my $clean_counts_arr_ref=sweep_low_counts_to_other(\@members_counts);
	compute_homogeneity_pvalue($clean_counts_arr_ref);

}

sub compute_homogeneity_pvalue{
	my $member_counts_arr_ref=shift;	# This is an array of an array references
	
	my $num_members=$#{$member_counts_arr_ref}+1;

	# Compute totals
	my @cat_total;
	my @member_totals_arr;
	my $all_total=0;
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		my @count_arr=@{${$member_counts_arr_ref}[$member_id]};
		print STDERR "Input: \n\t" . (join ",",@count_arr) . "\n";
		my $member_total=0;
		for(my $i=0; $i<=$#count_arr; $i++){
			$cat_total[$i]+=$count_arr[$i];
			$all_total+=$count_arr[$i];
			$member_total+=$count_arr[$i];
		}	
		push @member_totals_arr, $member_total;
	}
	print STDERR "Totals across all members:\n\t" . (join ",", @cat_total) . "\n";
	my $num_categories=$#cat_total+1;

	# Compute expected probabilities/percentages
	my @exp_perc;
	my @exp_perc_forshow;
	for(my $i=0; $i<$num_categories; $i++){
		$exp_perc[$i]+=$cat_total[$i]/$all_total;
		$exp_perc_forshow[$i]=sprintf("%2.4f", $exp_perc[$i]);
	}	
	print STDERR "Expected perc/probs across all members:\n\t" . (join ",", @exp_perc_forshow) . "\n";

	# Compute expected counts for each member
	print STDERR "Members totals:\n\t" . (join ",", @member_totals_arr) . "\n";
	my @members_expect_count_arr;
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		my @member_exp_count;
		for(my $i=0; $i<$num_categories; $i++){
			push @member_exp_count, int(.5+$exp_perc[$i]*$member_totals_arr[$member_id]);
		}
		print STDERR "Expected counts for member [$member_id]:\n\t" . (join ",", @member_exp_count) . "\n";
		push @members_expect_count_arr, \@member_exp_count;
	}

	# Merge 0 counts with lowest
	my ($members_merged_expected, $members_merged_observed)=
		merge_zero_expected_counts(\@members_expect_count_arr, $member_counts_arr_ref);

	# Compute the sum of the square differences between obs and exp
	my $chisqr_stat=0;
	my $zero_free_cat_counts;
	for(my $member_id=0; $member_id<$num_members; $member_id++){

		my @expected_counts=@{${$members_merged_expected}[$member_id]};
		my @observed_counts=@{${$members_merged_observed}[$member_id]};

		$zero_free_cat_counts=$#expected_counts+1;

		for(my $category_id=0; $category_id<$zero_free_cat_counts; $category_id++){
			my $diff=($observed_counts[$category_id]-$expected_counts[$category_id]);
			my $diffsqrd=$diff*$diff;
			my $chisq=$diffsqrd/$expected_counts[$category_id];
			$chisqr_stat+=$chisq;
		}
	}

	print STDERR "Chi-square stat: $chisqr_stat\n";
	my $df=($num_members-1)*($zero_free_cat_counts-1);
	print STDERR "Degrees of freedom = $df\n";
	my $uppertail=Statistics::Distributions::chisqrprob($df, $chisqr_stat);
	my $pvalue=1-$uppertail;
	print STDERR "P-value: $pvalue\n";
	print STDERR "if P-value < alpha, you accept that the members have the same proportions.\n";
	
	return($chisqr_stat, $df, $pvalue);
}

sub merge_zero_expected_counts{
	my $expected_counts_ref=shift;
	my $observed_counts_ref=shift;

	my $num_cat=$#{${$expected_counts_ref}[0]}+1;

	# Find where 0's are among all members
	my $num_members=$#{$expected_counts_ref}+1;
	my @zeros;
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		my @member_counts=@{${$expected_counts_ref}[$member_id]};
		for(my $cat=0; $cat<$num_cat; $cat++){
			if($member_counts[$cat]==0){
				$zeros[$cat]=1;
			}else{
				if(!defined($zeros[$cat])){
					$zeros[$cat]=0;
				}
			}
		}
	}
	print STDERR "Zeros: " . (join ",", @zeros) . "\n";

	# Sum up counts where they were 0's for any member
	my @merged_exp_cat;
	my @merged_obs_cat;
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		my @member_exp_counts=@{${$expected_counts_ref}[$member_id]};
		my @member_obs_counts=@{${$observed_counts_ref}[$member_id]};
		for(my $cat=0; $cat<$num_cat; $cat++){
			if($zeros[$cat]==1){
				$merged_exp_cat[$member_id]+=$member_exp_counts[$cat];
				$merged_obs_cat[$member_id]+=$member_obs_counts[$cat];
			}
		}
	}

	# Reconstruct obs/exp arrays without the zeros
	my @expected_arr;
	my @observed_arr;
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		my @member_counts=@{${$expected_counts_ref}[$member_id]};
		my @new_exp;
		my @new_obs;
		for(my $cat=0; $cat<$num_cat; $cat++){
			if($zeros[$cat]==0){
				push @new_exp, ${${$expected_counts_ref}[$member_id]}[$cat];
				push @new_obs, ${${$observed_counts_ref}[$member_id]}[$cat];
			}
		}	

		push @expected_arr, \@new_exp;
		push @observed_arr, \@new_obs;
	}

	# Find smallest category we can add the zero categories to
	my $new_length=$#{$expected_arr[0]}+1;
	my $min=${$expected_arr[0]}[0];
	my $min_idx=0;
	for(my $cat=0; $cat<$new_length; $cat++){
		if(${$expected_arr[0]}[$cat]<$min){
			$min=${$expected_arr[0]}[$cat];
			$min_idx=$cat;
		}	
	}	

	# Add counts to smallest category that has no zeros across all members
	for(my $member_id=0; $member_id<$num_members; $member_id++){
		${$expected_arr[$member_id]}[$min_idx]+=$merged_exp_cat[$member_id];
		${$observed_arr[$member_id]}[$min_idx]+=$merged_obs_cat[$member_id];

		#print STDERR "exp: " . (join ",", @{$expected_arr[$member_id]}) . "\n";
		#print STDERR "obs: " . (join ",", @{$observed_arr[$member_id]}) . "\n";

	}

	return(\@expected_arr, \@observed_arr);

}

sub sweep_low_counts_to_other{
	my $member_counts_arr=shift;

	# Determine which counts need to be swept together to remove 0's
	my @needs_sweeping;
	foreach my $count_arr_ref(@{$member_counts_arr}){
		my @count_arr=@{$count_arr_ref};
		for(my $i=0; $i<=$#count_arr; $i++){
			if($count_arr[$i]<1){
				$needs_sweeping[$i]++;
			}else{
				if(!defined($needs_sweeping[$i])){
					$needs_sweeping[$i]=0;
				}
			}
		}	
	}
	
	# Sweep up zero counts to last column if it doesn't sum to >0 
	my @all_new_counts_arr;
	my @other_arr;
	my $load_other=0;
	foreach my $count_arr_ref(@{$member_counts_arr}){

		my @count_arr=@{$count_arr_ref};
		my @new_count_arr;
		my $other=0;

		# Rebuild counts after remove those that need to be swept
		for(my $i=0; $i<=$#count_arr; $i++){
			if($needs_sweeping[$i]==0){
				push @new_count_arr, $count_arr[$i];
			}else{
				$other+=$count_arr[$i];
			}
		}	

		# keep track of whether there is anything and what to load to other.
		push @other_arr, $other;
		if($other>0){
			$load_other=1;
		}

		# Keep track of new zero free counts
		push @all_new_counts_arr, \@new_count_arr;
	}

	# If any of the other categories are 0, keep the other category
	if($load_other){
		for(my $member_id=0; $member_id<=$#all_new_counts_arr; $member_id++){
			push @{$all_new_counts_arr[$member_id]}, $other_arr[$member_id];
			print STDERR "Swept [$member_id]:\n\t" . (join ",", @{$all_new_counts_arr[$member_id]}) . "\n";
		}
	}

	return(\@all_new_counts_arr);
	
}

##########################################################################################

sub load_clusters{
	my $fn=shift;
	open(FH, "<$fn") || die "Could not open $fn\n";

	my %cluster_members;

	my $cur_clust_id="";
	while(<FH>){
		chomp;
		my ($clust_id, $member)=split /\t/, $_;
	
		if($clust_id ne ""){
			$cur_clust_id=$clust_id;
		}	

		if($cur_clust_id eq ""){
			die "Cluster ID could not be identified.\n";
		}
		push @{$cluster_members{$cur_clust_id}}, $member;

	}
	close(FH);

	foreach my $cluster_id(sort keys %cluster_members){
		print STDERR "$cluster_id\n";
		foreach my $member(@{$cluster_members{$cluster_id}}){
			print STDERR "\t$member\n";
		}
	}

	return(\%cluster_members);
}

##########################################################################################

sub load_summary_table{
	my $fn=shift;

	open(FH, "<$fn") || die "Could not open $fn\n";

	# read in the headers
	my $headers=<FH>;
	chomp $headers;

	my @classes=split /\t/, $headers;
	my %class_to_column_map;
	my %column_to_class_map;

	# remove the sample column
	shift @classes;

	my $i=0;
	print STDERR "Loaded count classes:\n";
	foreach my $class(@classes){
		$class_to_column_map{$class}=$i;
		$column_to_class_map{$i}=$class;
		$i++;
		print STDERR "\t$class -> $class_to_column_map{$class}\n";
	}
	
	# load the sample class counts
	my %sample_counts;
	while(<FH>){
		chomp;
		my @values=split /\t/, $_;
		my $sample_id=shift @values;
		my $total=shift @values;		# Toss these out
		$sample_counts{$sample_id}=\@values;
	}

	close(FH);
	
	return(\%column_to_class_map, \%class_to_column_map, \%sample_counts);	

}

##########################################################################################
