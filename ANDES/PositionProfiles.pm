#!/usr/bin/env perl

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################

package PositionProfiles;
use strict;
#use Math::CDF;

my @NUC_ARR =("A","C","G","T");
my @GAP_NUC_ARR=@NUC_ARR;
unshift @GAP_NUC_ARR, "-";

my @AA_ARR=("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");
my @GAP_AA_ARR=@AA_ARR;
unshift @GAP_AA_ARR, "-";


my $LN10=log(10);

use constant NUC_TYPE => "nuc";
use constant AA_TYPE => "aa";

sub new{
	my $this = shift;
        my $class = ref($this) || $this;
        my $self = {};
        bless $self, $class;
        $self->_initialize();
        return $self;
}

###############################################################################

sub get_alphabet_posprof{
	my $posprof_ref=shift;

	my @alphabet=sort keys %{$posprof_ref};
	#foreach my $char(@alphabet){print "$char\n";};
	
	# produce alphabet set without gap
	my @nogap_alphabet;
	foreach my $res (@alphabet){
		if($res ne "-" || $res ne ""){
			push @nogap_alphabet, $res;
		}
	}

	# Return complete and gapfree alphabet
	return(\@alphabet, \@nogap_alphabet);

}

sub get_alphabet_profarr{
	# Based on first element of position profile, extract the
	# alphabet
	my $profile_arr_ref=shift;

	my %alphabet_hash;

	# Pull alphabet
	for(my $i=0; $i<=$#{$profile_arr_ref}; $i++){
		foreach my $res (keys %{${$profile_arr_ref}[$i]}){
			$alphabet_hash{$res}=1;
		}
	}
	
	my @gap=sort keys %alphabet_hash;
	delete($alphabet_hash{"-"});
	my @nogap=sort keys %alphabet_hash;

	return(\@gap, \@nogap);
}

###############################################################################

sub add_gaps_to_profile{
# BOTH
# Given a profile array and a gapped sequence, a new profile array will be
# generated which contains gaps in it based on the where the gaps are in the 
# gapped sequence.  The beginning and end of the profile will have gaps set to 0.

	my $gapped_alignment_sequence_arr_ref=shift;
	my $gap_free_profile_arr_ref=shift;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($gap_free_profile_arr_ref);
	my @gap_res_arr=@{$gap_res_arr_ref};

	my @new_gapped_profile=();
	my @gapped_sequence_arr=@{$gapped_alignment_sequence_arr_ref};
	my $length=$#gapped_sequence_arr+1;
	my $gap_free_idx=0;

	for(my $i=0; $i<$length; $i++){

		my %prof;
		if($gapped_sequence_arr[$i] eq "-"){
			#print STDERR "Inserting gap at $i.\n";
			foreach my $nuc(@gap_res_arr){
				$prof{$nuc}=0;
			}
			$prof{"-"}=1;
		}else{
			foreach my $nuc(@gap_res_arr){
				$prof{$nuc}=${${$gap_free_profile_arr_ref}[$gap_free_idx]}{$nuc};
			}
			$gap_free_idx++;
		}

		$new_gapped_profile[$i]=\%prof;
	
	}

	# Remove gaps counts from ends
	for(my $i=0; $gapped_sequence_arr[$i] eq "-"; $i++){
		${$new_gapped_profile[$i]}{"-"}=0;
	}
	for(my $i=($length-1); $gapped_sequence_arr[$i] eq "-"; $i--){
		${$new_gapped_profile[$i]}{"-"}=0;
	}

	
	return(\@new_gapped_profile);
}

###############################################################################

sub compute_top_allele_frequencies{
# BOTH
	my $prof_arr_ref=shift;

	# Computes the primary and secondary most common alleles, and their frequencies
	# Returns array of alleles and frequences for 1st and 2nd most common allele.

	my $length=$#{$prof_arr_ref}+1;

	# For 1st most common allele
	my @alleles1;
	my @frequency1;
	# For 2nd most common allele
	my @alleles2;
	my @frequency2;

	for(my $i=0; $i<$length; $i++){
		my %freq_hash;
		my %pos_prof=%{${$prof_arr_ref}[$i]};
		
		# Bin nucleotides by frequency
		foreach my $nuc(@GAP_NUC_ARR){
			$freq_hash{$pos_prof{$nuc}}.=$nuc;
		}

		# Reverse sort by frequency
		my @sorted_freqs=sort {$b <=> $a} keys %freq_hash;

		# Save most common allele in array
		my $requested_freq=$sorted_freqs[0];
		push @alleles1, $freq_hash{$requested_freq};
		push @frequency1, $requested_freq;

		# Save 2nd most common allele in array
		my $requested_freq=$sorted_freqs[1];
		push @alleles2, $freq_hash{$requested_freq};
		push @frequency2, $requested_freq;
		
	}

	# Returns references to allele/frequency array
	return(\@alleles1, \@frequency1, \@alleles2, \@frequency2);

}

###############################################################################

sub normalized_shannon_entropy{
# BOTH
	my $prof_ref=shift;

	# Compute the shannon entropy contained at a given position
	my $num_nucs=0;
	my $sum=0;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($prof_ref);

	foreach my $nuc(@{$gap_res_arr_ref}){
		$sum+=${$prof_ref}{$nuc};
		$num_nucs++;
	}

	if($sum==0){
		print STDERR "Profile position has no information!\n";
		return (0);
	}

	my $sum_plog=0;
	foreach my $nuc(@{$gap_res_arr_ref}){
		my $p=${$prof_ref}{$nuc}/$sum;
		if($p>0){
			$sum_plog += $p*log($p);
		}
	}
	
	my $entropy=-$sum_plog;

	# Normalize so entropy is between 0 and 1
	my $normalization_factor=log($num_nucs);

	return($entropy/$normalization_factor);

}

###############################################################################

sub multiply_transitions{
# BOTH
	my $prof_src_ref=shift;
	my $prof_dst_ref=shift;

	# Computes the probability of transitioning from src profile to dst profile.
	# Pr(src)*Pr(dst)

	my %transition_hash;

	 my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($prof_src_ref);

	foreach my $src_nuc(@{$gap_res_arr_ref}){
		foreach my $dst_nuc(@{$gap_res_arr_ref}){

			if(${$prof_src_ref}{$src_nuc} eq "NaN" || ${$prof_dst_ref}{$dst_nuc} eq "NaN"){
				return(undef);
			}

			$transition_hash{$src_nuc}{$dst_nuc}=
				${$prof_src_ref}{$src_nuc} * ${$prof_dst_ref}{$dst_nuc};

		}
	}

	#foreach my $src_nuc(@GAP_NUC_ARR){
	#	foreach my $dst_nuc(@GAP_NUC_ARR){
	#		print STDERR $transition_hash{$src_nuc}{$dst_nuc} . "\t";
	#	}
	#	print STDERR "\n";
	#}
	#print STDERR "\n";

	return(\%transition_hash);
}

###############################################################################

sub compute_overall_transition_statistics{
# BOTH
	my $src_seq_id=shift;
	my $dst_seq_id=shift;
	my $gapped_profile_hash_ref=shift;

	#print STDERR "Comparing $a_seq_id vs. $b_seq_id\n";
	
	my $src_prof=${$gapped_profile_hash_ref}{$src_seq_id};
	my $dst_prof=${$gapped_profile_hash_ref}{$dst_seq_id};

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($src_prof);

	# Initialize cumulative transition hash
	my %cumulative_transition_hash;
	foreach my $src_nuc(@{$gap_res_arr_ref}){
		foreach my $dst_nuc(@{$gap_res_arr_ref}){
			$cumulative_transition_hash{$src_nuc}{$dst_nuc}=0.0;
		}
	}

	# Compute P(src & dst)
	my $prof_length=$#{$src_prof}+1;
	my $sum_score=0;

	my %src_freq=0;
	my %src_count=0;
	my $num_undefined_positions=0;

	for(my $i=0; $i<$prof_length; $i++){
		
		my $src_nuc_prof=normalize(${$src_prof}[$i]);
		my $dst_nuc_prof=normalize(${$dst_prof}[$i]);

		if(!defined($src_nuc_prof) || !defined($dst_nuc_prof)){
			print STDERR "Undefined position $i\n";
			$num_undefined_positions++;
			next;
		}
				
		my $transition_hash_ref=multiply_transitions($src_nuc_prof, $dst_nuc_prof);

		# Sum up the transition probabilities
		foreach my $src_nuc(@{$gap_res_arr_ref}){

			$src_freq{$src_nuc}+=${$src_nuc_prof}{$src_nuc};

			foreach my $dst_nuc(@{$gap_res_arr_ref}){
				$cumulative_transition_hash{$src_nuc}{$dst_nuc}+=
					${$transition_hash_ref}{$src_nuc}{$dst_nuc};
			}
		}
	}

	# Normalize to 0-1
	foreach my $src_nuc(@{$gap_res_arr_ref}){
		foreach my $dst_nuc(@{$gap_res_arr_ref}){	
			# Average the transitions across each position
			if($src_freq{$src_nuc} != 0 && defined($src_freq{$src_nuc})){
				$cumulative_transition_hash{$src_nuc}{$dst_nuc}/=$src_freq{$src_nuc};
			}else{
				$cumulative_transition_hash{$src_nuc}{$dst_nuc}=undef;
			}
		}
	}

	return(\%cumulative_transition_hash);
}

###############################################################################

sub rmse{
# BOTH
	# Compute root mean squared error
	my $posprof_a_ref=shift;
	my $posprof_b_ref=shift;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($posprof_a_ref);
	
	my $sum_squares=0;
	foreach my $nuc(@{$gap_res_arr_ref}){
		my $diff=(${$posprof_a_ref}{$nuc} - ${$posprof_b_ref}{$nuc});
		#print "my $diff=(${$prof_a_ref}{$nuc} - ${$prof_b_ref}{$nuc})\n";
		$sum_squares+=($diff*$diff);
	}
	return(sqrt($sum_squares/($#{$gap_res_arr_ref}+1)));
}

###############################################################################o

sub normalize{
	# Normalizes a single position
# BOTH
	my $posprof_ref=shift;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($posprof_ref);

	# Compute total count
	my $total=0;
	foreach my $nuc(@{$gap_res_arr_ref}){
		$total+=${$posprof_ref}{$nuc};
	}
	#print "total: $total\n";

	# normalize to between 0 and 1
	my %normalized_posprof;
	foreach my $nuc(@{$gap_res_arr_ref}){
		if($total>0){
			$normalized_posprof{$nuc}=${$posprof_ref}{$nuc}/$total;
		}else{
			return(undef);	
		}
	}

	return(\%normalized_posprof);
}

###############################################################################

sub normalize_gap_free{
	my $prof_ref=shift;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($prof_ref);	

	# Compute total count
	my $total=0;
	foreach my $res(@{$ungap_res_arr_ref}){
		$total+=${$prof_ref}{$res};
	}
	#print "total: $total\n";

	# normalize to between 0 and 1
	my %normalized_prof_hash;
	foreach my $res(@{$gap_res_arr_ref}){
		if($res eq "-"){
			$normalized_prof_hash{$res}=0;
		}else{
			if($total>0){
				$normalized_prof_hash{$res}=${$prof_ref}{$res}/$total;
				$normalized_prof_hash{$res}=sprintf("%6.5f", $normalized_prof_hash{$res});
			}else{
				$normalized_prof_hash{$res}="NaN";
			}
		}
	}

	return(\%normalized_prof_hash);
}

###############################################################################

sub normalize_profile{
	# Normalize the entire profile array
# BOTH
	my $prof_ref=shift;
	my @normalized_profile;

	for(my $i=0; $i<=$#{$prof_ref}; $i++){
		# If the normalized profile is undefined, we still want a place holder for it.
		$normalized_profile[$i]=normalize(${$prof_ref}[$i]);
	}

	return(\@normalized_profile);
}

###############################################################################

sub add_profile{
# BOTH
        my $new_prof_arr_ref=shift;             # array to add profile into
        my $aln_sequence_arr_ref=shift;         # how sequence aligned to consensus
        my $profile_arr_ref=shift;              # profile of sequence to add
        my $major_nuc_arr_ref=shift;            # sequence of profile to add

        my $new_idx;
        my $inp_idx=0;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($profile_arr_ref);

	my $prof_length=$#{$profile_arr_ref}+1;
        for($new_idx=0; $new_idx<=$#{$aln_sequence_arr_ref}; $new_idx++){


                # print "$new_idx: ${$aln_sequence_arr_ref}[$new_idx]\t${$major_nuc_arr_ref}[$inp_idx]\n";
                if(${$major_nuc_arr_ref}[$inp_idx] eq ${$aln_sequence_arr_ref}[$new_idx]){

			foreach my $nuc(@{$gap_res_arr_ref}){
				${${$new_prof_arr_ref}[$new_idx]}{$nuc}+=
					${${$profile_arr_ref}[$inp_idx]}{$nuc};
			}

                        $inp_idx++;
                }else{

			if($inp_idx>0 && $inp_idx<$prof_length){
				my $before_depth=sum_prof(${$profile_arr_ref}[$inp_idx-1]);
				my $after_depth=sum_prof(${$profile_arr_ref}[$inp_idx]);
				my $gap_depth=$before_depth>$after_depth?$after_depth:$before_depth;
				${${$new_prof_arr_ref}[$new_idx]}{"-"}+=$gap_depth;
			}
		}
        }
}

###############################################################################

sub compare_profiles{
# BOTH
	my $a_seq_id=shift;
	my $b_seq_id=shift;
	my $gapped_profile_hash_ref=shift;
	my $verbose_off=shift;
	
	my $warnings_off;
	if(defined($verbose_off)){
		$warnings_off=1;
	}else{
		$warnings_off=0;
	}

	#print STDERR "Comparing $a_seq_id vs. $b_seq_id\n";

	my $a_prof=${$gapped_profile_hash_ref}{$a_seq_id};
	my $b_prof=${$gapped_profile_hash_ref}{$b_seq_id};

	my $prof_length_A=$#{$a_prof}+1;
	my $prof_length_B=$#{$b_prof}+1;
	print STDERR "Profile Lengths A: $prof_length_A B: $prof_length_B\n";

	my $prof_length=$#{$a_prof}+1;
	my $sum_score=0;
	my $num_valid_values=0;
	for(my $i=0; $i<$prof_length; $i++){
		
		# Determine if positions are null, if either position is null,
		# then don't let it contribute to the RMSE calculation.
		my $a_null=0;
		my $b_null=0;
		if(!defined(${$a_prof}[$i])){
			$a_null=1;	
		}
		if(!defined(${$b_prof}[$i])){
			$b_null=1;
		}
		if($a_null || $b_null){
			if(!$warnings_off){
				if($a_null && $b_null){
					print STDERR "$i: $a_seq_id and $b_seq_id are NULL\n";
				}elsif($a_null){
					print STDERR "$i: $a_seq_id is NULL\n";
				}elsif($b_null){
					print STDERR "$i: $b_seq_id is NULL\n";
				}
			}
			next;
		}
		
		# Normalize
		my $a_nuc_prof=normalize(${$a_prof}[$i]);
		my $b_nuc_prof=normalize(${$b_prof}[$i]);

		if(!defined($a_nuc_prof) || !defined($b_nuc_prof)){
			next;
		}
				
		# Compute RMSE
		my $rmse_val=rmse($a_nuc_prof, $b_nuc_prof);

		$sum_score+=$rmse_val;
		$num_valid_values++;
	}

	my $prof_score=sprintf("%6.5f", ($sum_score/$num_valid_values));
	return($prof_score);
}

###############################################################################

sub compare_profiles_return_position_array{

	my $a_seq_id=shift;
	my $b_seq_id=shift;
	my $gapped_profile_hash_ref=shift;

	#print STDERR "Comparing $a_seq_id vs. $b_seq_id\n";
	my $a_prof=${$gapped_profile_hash_ref}{$a_seq_id};
	my $b_prof=${$gapped_profile_hash_ref}{$b_seq_id};

	my $rmse_arr_ref=compare_profiles_return_position_array_by_reference($a_prof, $b_prof);

	return($rmse_arr_ref);
}

###############################################################################

sub compare_profiles_return_position_array_by_reference{
# BOTH

	# Use RMSE to compare two profile arrays
	
	my $a_prof=shift;
	my $b_prof=shift;

	my $prof_length=$#{$a_prof}+1;
	my $sum_score=0;

	my @rmse_arr;
	for(my $i=0; $i<$prof_length; $i++){
		
		my $a_nuc_prof=normalize(${$a_prof}[$i]);
		my $b_nuc_prof=normalize(${$b_prof}[$i]);

		if(!defined($a_nuc_prof) || !defined($b_nuc_prof)){
			push @rmse_arr, undef;
		}else{
			my $rmse_val=rmse($a_nuc_prof, $b_nuc_prof);
			push @rmse_arr, $rmse_val;
		}
	}

	return(\@rmse_arr);
}

###############################################################################

sub readClustalAlnFile{
# Reads a clustal aln file into a hash of arrays(sequences), keyed by sequence id
	my $aln_file=shift;
	

	# Parse/Load the Aln file into a hash
	print STDERR "Reading ALN file...\n";
	
	open(ALN_FH, "<$aln_file") || die "Could not open $aln_file\n";
	
	my %seq_hash;
	my $clustal_file=0;
	while(<ALN_FH>){
		chomp;
		if($_=~/CLUSTAL/ || $_=~/MUSCLE/){
			$clustal_file=1;
			next;	
		}elsif($_=~/\*/){
			next;
		}else{
			my ($id, $sequence)=split /\s+/, $_;
			$sequence=uc($sequence);
			$sequence=~s/U/T/g;
			if($id ne ""){
				$seq_hash{$id}.=$sequence;
			}
		}
	}
	close(ALN_FH);
	
	if($clustal_file!=1){
		print STDERR "Are you sure this is clustal file?\n";
	}
	
	# Convert sequence into an array
	foreach my $id(keys %seq_hash){
		my @arr=split //, $seq_hash{$id};
		delete($seq_hash{$id});
		$seq_hash{$id}=\@arr;
	}

	print STDERR "done.\n";
	return(\%seq_hash);
	
}

################################################################################

sub load_prof{
	my $prof_name=shift;
	my @major_res_arr;
	my @profile_arr;

	open(FH, "<$prof_name") || die "Could not open profile $prof_name for reading.\n";

	my $alphabet_ok=0;
	my @alphabet;
	my $pos=0;
	while(<FH>){
		chomp;
		if(substr($_,0,1) eq "#"){
			# Look for magic number and load alphabet
			@alphabet=split /\t/, $_;
			my $magic_number=shift @alphabet;
			if($magic_number ne "#Major"){
				print STDERR  "Magic Number: $magic_number\n";
				die "Unexpected first line for defining profile alphabet.\n";
			}else{
				print STDERR "Identifiying profile alphabet type.\n";
				$alphabet_ok=1;
			}
		}else{
			if($alphabet_ok){
				my @line=split /\t/, $_;
				my $major_res=shift @line;

				if($#line > $#alphabet){
					die "More columns than alphabet supports.\n";
				}

				# Store each column in a hash according to alphabet
				for(my $i=0; $i<=$#line; $i++){
					${$profile_arr[$pos]}{$alphabet[$i]}=$line[$i];
				}
				# Store major residue
				push @major_res_arr, $major_res;
				$pos++;
			}else{
				die "Could not determine alphabet type.\n";
			}
		}
	}
	
	close(FH);

	return(\@major_res_arr, \@profile_arr);
}

################################################################################

sub compute_mode_residue_for_profile{
	my $prof_arr_ref=shift;

	my $prof_len=$#{$prof_arr_ref}+1;
	my $mode_sequence="";
	
	print STDERR "Profile Length: $prof_len\n";

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
	# Get the residues at this position
	my @res=@{$gap_res_arr_ref};

	# Finding major residue
	for(my $i=0; $i<$prof_len; $i++){

		# If the residue we initialize to is a gap, use the next one
		my $init_res=0;
		if($res[$init_res] eq "-"){
			$init_res++;
		}
		my $mode_var_res=$res[$init_res];
		my $mode_var_cnt=${${$prof_arr_ref}[$i]}{$res[$init_res]};

		# Go through each residue in the profile to see if there is another one with greater count
		foreach my $res(@res){
			if($mode_var_cnt < ${${$prof_arr_ref}[$i]}{$res} && $res ne "-"){
				$mode_var_cnt=${${$prof_arr_ref}[$i]}{$res};
				$mode_var_res=$res;
			}
		}

		# Appende the detected mode residue to the sequence
		$mode_sequence.=$mode_var_res;

	}

	return($mode_sequence);

}

################################################################################

sub output_prof{
# BOTH
	my $filename=shift;
	my $position_count_arr_ref=shift;
	my @major=split //, shift;
	my $residue_type=shift;

	my @pos_arr=@{$position_count_arr_ref};

	# Determine residue type
	my $alphabet;
	if($residue_type eq NUC_TYPE){
		$alphabet=\@GAP_NUC_ARR;	
	}elsif($residue_type eq AA_TYPE){
		$alphabet=\@GAP_AA_ARR;
	}else{
		my ($gap_alphabet, $nogap_alphabet)=get_alphabet_profarr($position_count_arr_ref);
		$alphabet=$gap_alphabet;
	}

	my $alphabet_len=$#{$alphabet}+1;

	# Output new output profile
	open(FH, ">$filename") || die "Could not open $filename for writing profile\n";

	# Output comment/header line
	my $alpha_str=join "\t", @{$alphabet};
	print FH "#Major\t$alpha_str\n";

	for(my $i=0; $i<=$#pos_arr; $i++){

		# Make all the undefined nucleotides 0
		foreach my $res(@{$alphabet}){
			if(!defined(${$pos_arr[$i]}{$res})){
				${$pos_arr[$i]}{$res}=0;
			}
		}

		# Concatenate counts together in order of the alphabet
		my @counts;
		for(my $letter_idx=0; $letter_idx<$alphabet_len; $letter_idx++){
			push @counts, ${$pos_arr[$i]}{${$alphabet}[$letter_idx]};
		}

		my $counts_str=join "\t", @counts;

		print FH "$major[$i]\t$counts_str\n";
	}

	close(FH);
}

################################################################################

sub input_fasta{
	my $filename=shift;
	open(FASTA, "<$filename") || die "Could not open $filename\n";

	my %fasta_hash;

	sub process_record{
		my $defline=shift;
		my $sequence=shift;
		my $hash_ref=shift;
		$defline=~s/^>//g;
		${$hash_ref}{$defline}=$sequence;	
	}

	my ($defline, $prev_defline, $sequence)=("","","");
	while(<FASTA>){
		chomp;

		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, $sequence, \%fasta_hash);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	process_record($prev_defline, $sequence, \%fasta_hash);
	close(FASTA);

	return(\%fasta_hash);

}

################################################################################

sub output_fasta{
        my $filename=shift; # output file name
	my $seqname=shift;  # sequence name to use in defline
        my $sequence=shift; # sequence
	my $width=shift;    # residues per line

        open(FH, ">$filename") || die "Could not open $filename\n";

        print FH ">$seqname\n";
        my $length=length($sequence);
        my $pos=0;
        do{
                my $out_width=($width>$length)?$length:$width;
                print FH substr($sequence, $pos, $width) . "\n";
                $pos+=$width;
                $length-=$width;
        }while($length>0);

        close(FH);
}

################################################################################

sub read_ClustalALN_into_profile{
	my $aln_file=shift;		# .aln file name from clustalw
	my $re_pattern=shift;		# regular expression pattern
	my %id_hash;
	my %all_id_hash;

	# returns a hash of ids, which each member contains the nucleotide profile
	
	print STDERR "Reading ALN file...\n";

	open(ALN_FH, "<$aln_file") || die "Could not open $aln_file\n";

	my $clustal_file=0;
	while(<ALN_FH>){
		chomp;
		if($_=~/CLUSTAL/){
			$clustal_file=1;
			next;
		}elsif($_=~/\*/){
			next;
		}else{
			my ($id, $sequence)=split /\s+/, $_;
			if($id ne ""){
				$sequence=uc($sequence);
				$all_id_hash{$id}=1;
				if($id=~/$re_pattern/){
					$id_hash{$id}.=$sequence;
				}
			}
		}
	}
	close(ALN_FH);

	if($clustal_file!=1){
		print STDERR "Are you sure this is clustal file?\n";
	}

	print STDERR "done.\n";

	my $num_kept=keys %id_hash;
	my $num_ids=keys %all_id_hash;
	print STDERR "Pattern: /$re_pattern/ $num_kept/$num_ids\n";

	# Generate profile based on sequences
	my $profile_ref=sequence_hash_into_profile(\%id_hash);

	return($profile_ref);	
}

################################################################################

sub sequence_hash_into_profile{
	my $seq_hash_ref=shift;
	my $residue_type=shift;

	print STDERR "Converting aligned/gapped sequences into profile...\n";

	my @profile;	# this is returned
	my $length=-1;

	foreach my $id(keys %{$seq_hash_ref}){

		print STDERR ">$id\n";
		my @seq_arr=@{${$seq_hash_ref}{$id}};
		
		#print STDERR (join "", @seq_arr) . "<EOS>\n";
		# Remove leading and trailing -, because they are not really gaps.
		# They are just incomplete sequences.

		my $seq_len=$#seq_arr+1;
		for(my $i=0; $i<$seq_len; $i++){
			if($seq_arr[$i] eq "-"){
				$seq_arr[$i]=" ";
			}else{
				last;
			}
		}

		for(my $i=($seq_len-1); $i>=0; $i--){
			if($seq_arr[$i] eq "-"){
				$seq_arr[$i]=" ";
			}else{
				last;
			}
		}
		#print STDERR (join "", @seq_arr) . "<EOS>\n";

		# Assume length of consensus is the length of the first sequence, but die if not all
		#  lengths are the same length.
		my $curlength=$#seq_arr+1;
		if($length==-1){
			$length=$curlength;
			for(my $i=0; $i<$length; $i++){
				$profile[$i]=init_pos_prof($residue_type);
			}
		}elsif($length != $curlength){
			die "Error, for some reason the length of alignment for $id is not the same as all the other sequences ($length != $curlength)\n";
		}

		# Compute the profile at each position
		#       Each position is represented in an array
		#       Each nucleotide code is summed up in a hash
		for(my $i=0; $i<$length; $i++){
			if($seq_arr[$i] ne " "){
				ambig_to_prof($seq_arr[$i], \%{$profile[$i]}, $residue_type);
			}
		}
	}

	print STDERR "\ndone.\n";
	return(\@profile);

}

################################################################################

sub rebuild_filtered_all_gaps{
	my @profile_references=@_;
	my $num_profiles=$#profile_references+1;

	my @new_profile_references=();

	# This function will go through all the aligned profiles and then
	# remove positions where each position is 100% gaps across all
	# profiles.
	
	# Get the length for all profiles
	my $prof_length=$#{$profile_references[0]};
	foreach my $prof_ref(@profile_references){
		if($#{$profile_references[0]}!=$prof_length){
			print STDERR "$#{$profile_references[0]}}!=$prof_length\n";
			die "Profile lengths are not all the same.\n";
		}
	}
	$prof_length++;


	# For each position determine if there are any gap only positions
	my $new_prof_index=0;
	for(my $i=0; $i<$prof_length; $i++){

		# Look across all profiles at specific position for gap only's
		my $nuc_found=0;
		for(my $prof_index=0; $prof_index<$num_profiles; $prof_index++){
			my $pos_prof_ref=${$profile_references[$prof_index]}[$i];
			if(!gap_only($pos_prof_ref)){
				$nuc_found++;	
			}
		}

		# If gap only, then exclude it from filtered profiles.
		if($nuc_found){
			for(my $prof_index=0; $prof_index<$num_profiles; $prof_index++){
				${$new_profile_references[$prof_index]}[$new_prof_index]=${$profile_references[$prof_index]}[$i];
			}
			$new_prof_index++;
		}
	}

	print STDERR "Incoming Length: $prof_length  All Gap Filtered Length: $new_prof_index\n";

	return(@new_profile_references);

}

################################################################################

sub filter_prof_depth{
	my $prof_arr_ref=shift;
	my $min_depth=shift;
	my @new_profile=();

	# If the sum of the depth is below the min depth, zero out all the residue counts

	print STDERR "Filtering profile by min depth ($min_depth)...\n";

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);

	my $prof_length=$#{$prof_arr_ref}+1;

	for(my $i=0; $i<$prof_length; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];

		my $tot_depth=0;
		foreach my $res(@{$gap_res_arr_ref}){
			$tot_depth += ${$pos_prof_ref}{$res};
		}	

		if($tot_depth<$min_depth){
			foreach my $res(@{$gap_res_arr_ref}){
				${$new_profile[$i]}{$res}=0;
			}
		}else{
			foreach my $res(@{$gap_res_arr_ref}){
				${$new_profile[$i]}{$res}=${$pos_prof_ref}{$res};
			}
		}
	}

	return(\@new_profile);
}

################################################################################

sub filter_prof_basecount{
	my $prof_arr_ref=shift;
	my $base_count_thres=shift;
	my @new_profile=();

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);

	my $prof_length=$#{$prof_arr_ref}+1;

	for(my $i=0; $i<$prof_length; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];
		foreach my $res(@{$gap_res_arr_ref}){
			if(${$pos_prof_ref}{$res}<$base_count_thres){
				${$new_profile[$i]}{$res}=0;
			}else{
				${$new_profile[$i]}{$res}=${$pos_prof_ref}{$res};
			}
		}	
	}

	return(\@new_profile);
}

################################################################################

sub filter_prof_percentage{
	my $prof_arr_ref=shift;
	my $base_perc_thres=shift;
	my @new_profile=();

	my $prof_length=$#{$prof_arr_ref}+1;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);

	for(my $i=0; $i<$prof_length; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];
		
		# Count up the total counts
		my $pos_count=0;
		foreach my $res(@{$gap_res_arr_ref}){
			$pos_count+=${$pos_prof_ref}{$res};
		}


		# Go through each nucleotide count
		foreach my $res(@{$gap_res_arr_ref}){

			my $res_perc;

			if($pos_count>0){
				$res_perc=${$pos_prof_ref}{$res}/$pos_count;
			}else{
				$res_perc=0;
			}

			if($res_perc<$base_perc_thres){
				${$new_profile[$i]}{$res}=0;
			}else{
				${$new_profile[$i]}{$res}=${$pos_prof_ref}{$res};
			}
		}	

	}

	return(\@new_profile);
}

################################################################################

#sub filter_prof_binomial{
#
#	# For each nucleotide in the distribution, removes the allele if the probability
#	# of that nucleotide being greater than 0 is better than the threshold
#	# specified, if resampled at the same depth were to be performed.
#
#	# In other words, if we were to resample the same position again, if the probability
#	# of that allele being nonzero is greater than the specified threshold, keep it.
#
#	my $prof_arr_ref=shift;
#	my $conf_thres=shift; # The higher the number, the greater the stringency.
#	my @new_profile=();
#
#	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
#
#	my $prof_length=$#{$prof_arr_ref}+1;
#
#	for(my $i=0; $i<$prof_length; $i++){
#		my $pos_prof_ref=${$prof_arr_ref}[$i];
#		
#		# Count up the total counts
#		my $pos_count=0;
#		foreach my $res(@{$gap_res_arr_ref}){
#			$pos_count+=${$pos_prof_ref}{$res};
#		}
#
#		# Compute prob of each nucleotide
#		my %res_perc;
#		foreach my $res(@{$gap_res_arr_ref}){
#			if($pos_count>0){
#				$res_perc{$res}=${$pos_prof_ref}{$res}/$pos_count;
#			}else{
#				$res_perc{$res}=0;
#			}
#		}
#
#		# Test each nucleotide for inclusion
#		my $new_total=0;
#		foreach my $res(@{$gap_res_arr_ref}){
#
#			my $p=$res_perc{$res};
#			#print STDERR "Count = $pos_count p($nuc)=$p\n";
#			my $prob_not_zero=1-Math::CDF::pbinom(0, $pos_count, $p);
#			#print STDERR "P(X>0)=$prob_not_zero\n";
#
#			if($prob_not_zero>$conf_thres){
#				${$new_profile[$i]}{$res}=${$pos_prof_ref}{$res};
#				$new_total+=${$pos_prof_ref}{$res};
#			}else{
#				${$new_profile[$i]}{$res}=0;
#			}
#		}	
#
#		# If all counts were filtered, revert and report.
#		if($new_total==0){
#			#print STDERR "Warning: Position $i could not be filtered.  All counts were insignificant.\n";
#			foreach my $res(@{$gap_res_arr_ref}){
#				${$new_profile[$i]}{$res}=${$pos_prof_ref}{$res};
#			}
#		}
#	}
#
#	return(\@new_profile);
#}

################################################################################

sub gap_only{
	my $prof_ref=shift;
	
	# Returns 1 if the position has no bases in it

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($prof_ref);

	my $sum=0;
	foreach my $res(@{$gap_res_arr_ref}){
		if($res ne "-"){
			$sum+=${$prof_ref}{$res};
		}
	}

	return($sum==0);
}

################################################################################

sub freq_at_pos{
	my $prof_arr_ref=shift;
	my $pos=shift;
	my $res=shift;
	my $normalized=shift;

	# For the specified res and position, it will grab the frequency
	# or the probability of that res in that position.  If you specify true
	# for normalized, the returned value will be the probability, otherwise
	# it will be the frequency (count).

	my $pos_prof_ref=${$prof_arr_ref}[$pos];
	if($normalized){
		$pos_prof_ref=normalize($pos_prof_ref);
		if(!defined($pos_prof_ref)){
			return("NaN");
		}
	}

	return(${$pos_prof_ref}{uc($res)});
}

################################################################################

sub find_residue_pos_above_threshold{
	my $prof_arr_ref=shift;
	my $threshold=shift;
	my $residue=shift;
	
	# This function will report an array of all the positions that contain
	# the specified nucleotide equaling or exceeding the specified threshold.

	my @target_res_pos_arr;

	my $prof_len=$#{$prof_arr_ref}+1;
	for(my $i=0; $i<$prof_len; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];
		my $norm_pp=normalize($pos_prof_ref);
		if(defined($norm_pp)){
			if(${$norm_pp}{$residue} >= $threshold){
				push @target_res_pos_arr, $i;
			}
		}
	}

	return(\@target_res_pos_arr);

}

################################################################################

sub filter_gaps{
	my $prof_arr_ref=shift;
	my $threshold=shift;

	# Takes the input profile array and filters out positions that contain
	# gaps greater than the specified threshold.

	my @new_prof_arr;

	my $prof_len=$#{$prof_arr_ref}+1;

	for(my $i=0; $i<$prof_len; $i++){
		my $gap_prob=freq_at_pos($prof_arr_ref, $i, "-");
		if($gap_prob <= $threshold){
			push @new_prof_arr, ${$prof_arr_ref}[$i];
		}
	}

	return(\@new_prof_arr);

}

################################################################################

sub remove_pos_with_all_gaps{
	my $prof_arr_ref=shift;

	# If the position is 100% gaps, remove it
	
	my @new_prof_arr;
	my $prof_len=$#{$prof_arr_ref}+1;
	
	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
	
	for(my $i=0; $i<$prof_len; $i++){

		my $res_counts=0;
		foreach my $res(@{$ungap_res_arr_ref}){
			$res_counts+=${${$prof_arr_ref}[$i]}{$res};
		}

		my $gap_count=${${$prof_arr_ref}[$i]}{"-"};

		if($res_counts>0 || $gap_count==0){
			push @new_prof_arr, ${$prof_arr_ref}[$i];	
		}
	}
	return(\@new_prof_arr);
}


################################################################################

sub estimate_qv{
	my $prof_ref=shift;
	my $round=shift;

	# Estimates the quality value for a nuc profile, assuming that the most
	# common allele is the correct one.

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($prof_ref);

	my $prob_major_allele=0;
	my $major_allele="";

	my $norm_pos_prof=normalize($prof_ref);
	foreach my $res(@{$gap_res_arr_ref}){
		if(${$norm_pos_prof}{$res}>$prob_major_allele){
			$prob_major_allele=${$norm_pos_prof}{$res};
			$major_allele=$res;
		}	
	}

	my $error_rate=1-$prob_major_allele;

	my $qv;
	if($error_rate==0){
		$qv=40;
	}else{
		$qv=-10*log($error_rate)/$LN10;
		if($round){
			$qv=int($qv+0.5);
		}else{
			$qv=sprintf("%5.4f",$qv);
		}
	}

	return($qv);

}

################################################################################

sub freq_at_pos{
	my $prof_arr_ref=shift;
	my $pos=shift;
	my $nuc=shift;

	my $normalized_pos_prof=normalize(${$prof_arr_ref}[$pos]);
	return(${$normalized_pos_prof}{$nuc});
	
}

################################################################################

sub compute_residue_composition{
	my $prof_arr_ref=shift;
	my $begin=shift;
	my $end=shift;
	my $max_gap_threshold=shift;

	# This function will compute the nucleotide composition of the profile between
	# begin and end.  Returns a hash with the normalized counts.  If
	# the gap threshold is exceed for a certain position, that position
	# will be ignored, otherwise, the position will be renormalized without the
	# gap, so the nucleotides sum up to 1.0.  The idea is that if the gap representation
	# is too high, then that position is probably a sequencing error.

	my $prof_len=$#{$prof_arr_ref}+1;
	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
	
	# Make sure we get what we want
	if($begin>$end){
		($begin,$end)=($end,$begin);
	}
	$begin=($begin<0)? 0:$begin;
	$end=($end>$prof_len)? $prof_len:$end;
	if($max_gap_threshold>1){
		die "\$max_gap_threshold for compute_residue_content(), must be between 0-1.\n";
	}

	# Initialize the hash
	my %sum_hash;
	foreach my $nuc(@{$gap_res_arr_ref}){
		$sum_hash{$nuc}=0;
	}
	my $valid_positions=0;

	# Sum up the nucs and valid positions
	for(my $i=$begin; $i<$end; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];
		my $norm_pp_ref=normalize($pos_prof_ref);
		if(defined($norm_pp_ref)){
			# Check the prob of gap.  Make sure it is less than the threshold.
			if(${$norm_pp_ref}{"-"} <= $max_gap_threshold){
				$valid_positions++;
				my $gap_free_ref=normalize_gap_free($norm_pp_ref);			
				foreach my $res(@{$gap_res_arr_ref}){
					$sum_hash{$res}+=${$gap_free_ref}{$res};
				}
			}
		}
	}
		
	# Compute the average for the window
	foreach my $res(@{$gap_res_arr_ref}){
		$sum_hash{$res}/=$valid_positions;
	}

	return(\%sum_hash);

}

################################################################################

sub compute_total_residues_assayed{
	my $prof_arr_ref=shift;
	my $begin=shift;
	my $end=shift;

	# This function will compute the total number of bases sequenced between
	# begin and end.  It will return the overall and per base.

	my $prof_len=$#{$prof_arr_ref}+1;
	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
	
	# Make sure we get what we want
	if($begin>$end){
		($begin,$end)=($end,$begin);
	}
	$begin=($begin<0)? 0:$begin;
	$end=($end>$prof_len)? $prof_len:$end;

	# Initialize the hash/variables
	my %sum_hash;
	foreach my $res(@{$gap_res_arr_ref}){
		$sum_hash{$res}=0;
	}
	my $total=0;

	# Sum up the nucs and valid positions within the window
	for(my $i=$begin; $i<$end; $i++){
		my $pos_prof_ref=${$prof_arr_ref}[$i];
		foreach my $res(@{$gap_res_arr_ref}){
			$sum_hash{$res}+=${$pos_prof_ref}{$res};
			$total+=${$pos_prof_ref}{$res};
		}
	}

	return($total, \%sum_hash);
}

################################################################################

sub is_indel_from_homopolymer{
	my $pos=shift;			# Which position to analyze
	my $base_of_interest=shift;	# Which nucleotide to test for up/down stream track
	my $prof_arr_ref=shift;		# Reference to profile array
	my $min_threshold=shift;	# Min threshold for determining end of track (between 0-1)

	# Given a position in the profile array, will check to see if the base
	# of interest is part of an up or downstream homopolymer track.

	my @dir_arr=(-1, 1);	# Check 5' and 3' of position
	my $track_len=0;
	my $prof_arr_len=$#{$prof_arr_ref}+1;

	foreach my $dir (@dir_arr){
		my $end_of_track=0;
		my $cur_pos=$pos+$dir;
		while($cur_pos>=0 && $cur_pos<$prof_arr_len && !$end_of_track){
			my $pos_freq=freq_at_pos($prof_arr_ref, $cur_pos, $base_of_interest, 1);
			if($pos_freq>=$min_threshold && $pos_freq!=0){
				$track_len++;
			}else{
				$end_of_track=1;
			}	
			$cur_pos+=$dir;
		}
	}

	return($track_len);
}

################################################################################

sub prof_to_ambig{
	my $pos_prof_ref=shift;
	my $nuc;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($pos_prof_ref);
	my $alphabet_length=$#{$ungap_res_arr_ref}+1;
	if($alphabet_length>4){
		print STDERR "Alphabet length is: $alphabet_length residues.\n";
		die "Profile appears not to be of nucleic acids.\n";
	}

	my %nuc_hash;
	foreach my $nuc(@GAP_NUC_ARR){
		if($nuc eq "-"){
			next;
		}
		my $nuc_count=${$pos_prof_ref}{$nuc};
		if(defined($nuc_count) && $nuc_count>0){
			$nuc_hash{$nuc}=1;
		}
	}
	my $code=get_amb(\%nuc_hash);

	return($code);
}

#------------------------------------------------------------------------------

sub get_amb{
        my $nuc_hash_ref=shift;
        my $code;

        # The order of these if statements are important.

        if(has($nuc_hash_ref, "A")){
                $code="A";
        }
        if(has($nuc_hash_ref, "C")){
                $code="C";
        }
        if(has($nuc_hash_ref, "G")){
                $code="G";
        }
        if(has($nuc_hash_ref, "T")){
                $code="T";
        }
        if(has($nuc_hash_ref, "AC")){
                $code="M";
        }
        if(has($nuc_hash_ref, "AG")){
                $code="R";
        }
        if(has($nuc_hash_ref, "AT")){
                $code="W";
        }
        if(has($nuc_hash_ref, "CG")){
                $code="S";
        }
        if(has($nuc_hash_ref, "CT")){
                $code="Y";
        }
        if(has($nuc_hash_ref, "GT")){
                $code="K";
        }
        if(has($nuc_hash_ref, "ACG")){
                $code="V";
        }
        if(has($nuc_hash_ref, "ACT")){
                $code="H";
        }
        if(has($nuc_hash_ref, "AGT")){
                $code="D";
        }
        if(has($nuc_hash_ref, "CGT")){
                $code="B";
        }
        if(has($nuc_hash_ref, "GATC")){
                $code="N";
        }
	if(!has($nuc_hash_ref, "A") &&
	   !has($nuc_hash_ref, "T") &&
	   !has($nuc_hash_ref, "G") &&
	   !has($nuc_hash_ref, "C")){
		$code="N";
		print STDERR "Warning: Position has no nucleotides defined. Assuming N.\n";
	}
        return($code);
}

sub has{
        # See if all the members in the string, have a representive in the hash
        my $hash_ref=shift;
        my %hash=%{$hash_ref};  # Hash of nucleotide counts
        my $str=shift;          # String containing representatives of ambiguity code

        # Test all the nucs in the hash to see if they exist in the hash
        my @nucs=split //, $str;
        my $num_hits=0;
        foreach my $nuc(@nucs){
                if(defined($hash{$nuc})){
                        $num_hits++;
                }
        }

        # Only return true if all the members in the string are represented in the hash
        if($num_hits==($#nucs+1)){
                return(1);
        }else{
                return(0);
        }
}

################################################################################

sub ambig_to_prof_aa{
	my $aa=shift;
	my $prof_ref=shift;

	if(
		$aa eq "B"
	){
		${$prof_ref}{"D"}+=.5;
		${$prof_ref}{"N"}+=.5;
	}elsif(
		$aa eq "Z"
	){
		${$prof_ref}{"E"}+=.5;
		${$prof_ref}{"Q"}+=.5;
	}elsif(
		$aa eq "J"
	){
		${$prof_ref}{"I"}+=.5;
		${$prof_ref}{"L"}+=.5;
	}elsif(
		# Invalid amino acids
		!($aa eq "O" ||			# O possibly Pyrolysine
		$aa eq "U")			# U possibly Selenocysteine
	){

		if($aa ne "X"){ # X is valid, but we won't count it
			${$prof_ref}{$aa}+=1;
		}
	}else{
		die "Invalid AA: ($aa)\n";
	}
}

sub ambig_to_prof_nuc{
	
	my $nuc=shift;
	my $prof_ref=shift;

	#--------------------------------------
	# Nonambiguous

	if(
		$nuc eq "A" ||
		$nuc eq "T" ||
		$nuc eq "G" ||
		$nuc eq "C" ||
		$nuc eq "-"
	){
		${$prof_ref}{$nuc}+=1;
	}

	#--------------------------------------
	# 2-way ambiguous

	if(
		$nuc eq "M" ||
		$nuc eq "R" ||
		$nuc eq "W"
	){
		${$prof_ref}{"A"}+=0.5;
	}
	if(
		$nuc eq "W" ||
		$nuc eq "Y" ||
		$nuc eq "K"
	){
		${$prof_ref}{"T"}+=0.5;
	}
	if(
		$nuc eq "R" ||
		$nuc eq "S" ||
		$nuc eq "K"
	){
		${$prof_ref}{"G"}+=0.5;
	}
	if(
		$nuc eq "M" ||
		$nuc eq "S" ||
		$nuc eq "Y"
	){
		${$prof_ref}{"C"}+=0.5;
	}

	#--------------------------------------
	# 3-way ambiguous

	if(
		$nuc eq "V" ||
		$nuc eq "H" ||
		$nuc eq "D"
	){
		${$prof_ref}{"A"}+=0.3333;
	}
	if(
		$nuc eq "H" ||
		$nuc eq "D" ||
		$nuc eq "B"
	){
		${$prof_ref}{"T"}+=0.3333;
	}
	if(
		$nuc eq "V" ||
		$nuc eq "D" ||
		$nuc eq "B"
	){
		${$prof_ref}{"G"}+=0.3333;
	}
	if(
		$nuc eq "V" ||
		$nuc eq "H" ||
		$nuc eq "B"
	){
		${$prof_ref}{"C"}+=0.3333;
	}

	if($nuc=~/[EFIJLOPQZ]/){
		die "Invalid nucleic acid: $nuc\n";
	}

	#--------------------------------------
	# N's and X's are will not contribute any information.

}

################################################################################

sub init_pos_prof{
	my $res_type=shift;
	my %prof_hash;

	my $res;
	if($res_type eq AA_TYPE){
		foreach my $res(@GAP_AA_ARR){
			$prof_hash{$res}=0;
		}	
	}elsif($res_type eq NUC_TYPE){
		foreach my $res(@GAP_NUC_ARR){
			$prof_hash{$res}=0;
		}	
	}
	return(\%prof_hash);
}

################################################################################

sub ambig_to_prof{
	# This makes calls to ambig_to_prof_nuc and ambig_to_prof_aa

	my $res=shift;
	my $prof_ref=shift;
	my $res_type=shift;

	if($res_type eq AA_TYPE){
		ambig_to_prof_aa($res, $prof_ref);
	}elsif($res_type eq NUC_TYPE){
		ambig_to_prof_nuc($res, $prof_ref);
	}else{
		die "Error: In call to ambig_to_prof, residue type not specified. ($res_type)\n";
	}

}

################################################################################

sub print_pos_prof{
	my $pos_prof_ref=shift;
	my $fh=shift;

	if(!defined($fh)){
		$fh=\*STDOUT;
	}

	#print STDERR join ",",(keys %{$pos_prof_ref});

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($pos_prof_ref);

	foreach my $res(@{$gap_res_arr_ref}){
		my $val_str=sprintf("%4.4f", ${$pos_prof_ref}{$res});
		print {$fh} "$res: $val_str\t";
	}

	# Output place holder for null position
	if($#{$gap_res_arr_ref}==-1){
		print {$fh} "[NULL]: ";
	}

	print {$fh} "\n";
}

################################################################################

sub print_prof{
	my $position_count_arr_ref=shift;
	my $fh=shift;
	my $length=$#{$position_count_arr_ref}+1;

	for(my $i=0; $i<$length; $i++){
		print_pos_prof(${$position_count_arr_ref}[$i], $fh);
	}
}

################################################################################

sub sum_prof{

	# Computes the depth of evidence at the position.  ie. includes gaps
	my $pos_prof_ref=shift;
	my $depth=0;

	#my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_posprof($pos_prof_ref);

	foreach my $nuc(keys %{$pos_prof_ref}){
		$depth+=${$pos_prof_ref}{$nuc};
        }
	return($depth);
}

################################################################################

sub count_null_positions{
	my $prof_arr_ref=shift;
	my ($num_null, $num_valid)=(0,0);

	foreach my $pos_ref(@{$prof_arr_ref}){
		my $depth=sum_prof($pos_ref);
		if($depth>0){
			$num_valid++;
		}else{
			$num_null++;
		}
	}

	return($num_null, $num_valid);

}

################################################################################

sub profile_arr_length{
	my $prof_arr_ref=shift;
	my $prof_length=$#{$prof_arr_ref}+1;
	return($prof_length);
}

################################################################################

sub count_ambiguous_positions_across_profile{
        my $prof_arr_ref=shift;
	my $ambiguous_count=0;

        my $prof_length=$#{$prof_arr_ref}+1;

	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);

        for(my $i=0; $i<$prof_length; $i++){

		my $non_zero_count=0;

		my $prof_ref=${$prof_arr_ref}[$i];
		
                foreach my $res(@{$gap_res_arr_ref}){
			if(${$prof_ref}{$res}>0){
				$non_zero_count++;
			}
                }

		if($non_zero_count>1){
			$ambiguous_count++;
		}

        }

        return($ambiguous_count);
}

################################################################################

sub is_null_prof{
	my $pos_prof_ref=shift;
	my $sum=sum_prof($pos_prof_ref);
	if($sum==0){
		return(1);
	}else{
		return(0);
	}
}

################################################################################

sub size_null_ends{
	my $prof_arr_ref=shift;
	my $prof_len=profile_arr_length($prof_arr_ref);

	# Look for null from beginning of array
	my $begin_nulls=0;
	for(my $i=0; $i<$prof_len; $i++){
		if(is_null_prof(${$prof_arr_ref}[$i])){
			$begin_nulls++;
		}else{
			last;
		}
	}

	# Look for null from end of array
	my $end_nulls=0;
	for(my $i=($prof_len-1); $i>=0; $i--){
		if(is_null_prof(${$prof_arr_ref}[$i])){
			$end_nulls++;
		}else{
			last;
		}
	}
	
	return($begin_nulls, $end_nulls);

}

################################################################################

sub trim{
	my $prof_arr_ref=shift;
	my $start=shift;	# Keep start from 0
	my $end=shift;		# Keep end

	my $prof_len=profile_arr_length($prof_arr_ref);
	my @new_prof_array;
	
	my $new_arr_idx=0;
	for(my $i=$start; $i<=$end; $i++){
		$new_prof_array[$new_arr_idx]=${$prof_arr_ref}[$i];
		$new_arr_idx++;
	}
	
	return(\@new_prof_array);
}

################################################################################

sub mask_positions_over_threshold{
	# The prof_arr should already be normalized
	my $prof_arr_ref=shift;
	my $threshold=shift;

	print STDERR "Thres: $threshold\n";
	
	my $prof_len=$#{$prof_arr_ref}+1;
 	my ($gap_res_arr_ref, $ungap_res_arr_ref)=get_alphabet_profarr($prof_arr_ref);
	my @threshold_mask;

	for(my $i=0; $i<$prof_len; $i++){

                my $thres_exceeded=0;
                my $prof_ref=${$prof_arr_ref}[$i];

                foreach my $res(@{$gap_res_arr_ref}){
		#	print STDERR "$res: ${$prof_ref}{$res} ";
                        if(${$prof_ref}{$res} >= $threshold){
				$thres_exceeded=1;
			}
                }

		#print STDERR " --> $thres_exceeded\n";
	
		push @threshold_mask, $thres_exceeded;

	}

	return(\@threshold_mask);

}

################################################################################

1;
