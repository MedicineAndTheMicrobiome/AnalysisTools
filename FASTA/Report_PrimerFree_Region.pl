#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_e $opt_E $opt_f $opt_p $opt_m $opt_o $opt_a);
getopts("e:E:f:p:m:o:a:");

# Need at least EMBOSS 5.0, 4 is broken
my $FUZZ_NUC_BIN="/usr/local/packages/EMBOSS-5.0.0/bin/fuzznuc";
my $fuzznuc_found_str=`which $FUZZ_NUC_BIN`;
if($fuzznuc_found_str=~/Command not found/){
	print STDERR "$FUZZ_NUC_BIN not found.\n";
}

# Defaults
my $ERROR=6;
my $MARGIN=100;

my $usage="
	
	$0 

	-f <reads multifasta file>
	-p <primer multifasta file>
	-o <output root name>

	[-e <acceptable maximum error rate (excludes ambiguity codes), def. $ERROR>]

	[-E <acceptable percentage identity (excludes ambiguity codes), ie. 97>]
		(Note: this won't be exact if your primers are too short)

	Only one of -e or -E may be specified.

	[-m <margin from ends to allow primers to be believed, def $MARGIN bps>]
		(This should be the distance from the M13 primers)

	[-a <estimated amplicon length>]
		(If you use this option, then you probably only want to have one pair
		of primers per primer multi-fasta file, or else, it will assume that 
		all your amplicons are the same length.  By using this option you will
		increase the amount of margin given to the 3' end of the read dynamically
		based on the read length.)
	
	
	Reads in a multi-fasta of reads and a multifasta of primers.

	Returns a list of trims in <output_root>.clear:

		<read_id>\\t<begin>\\t<end>\\t<primername>/<mismatches>/<five or three prime>[;]<primername>/<mismatches>/<five or three prime>\\n
		...

	Also produces a statistics file, <output_root>.stats:
		<filter.fieldname>=<value>
		<filter.fieldname>=<value>
		...
		<filter.fieldname>=<value>


	This program uses fuzznuc: $FUZZ_NUC_BIN

	You should use this output to trim the fasta and quality sequences.


	------M13For----------16sFor=========================16sRev------------M13Rev------------
	
	      |------------------------Read------------------------------------------------|

	                            |=========Amplicon======|

	      |-----------Margin--------|                |----Margin + (ReadLen-AmpLen)----|
	
";


if(
	!defined($opt_f) ||
	!defined($opt_p) ||
	!defined($opt_o)
){
	die $usage;
}

if(defined($opt_e) && defined($opt_E)){
	die "You can't define both error rate and error percentage.\n";
}

my $error_rate=(!defined($opt_e)?$ERROR:$opt_e);
my $min_perc_id=$opt_E/100.0;
if(defined($opt_E) && $min_perc_id>100.0){
	die "The percent identity that you specified exceeds 100%.\n";
}

my $reads_fasta=$opt_f;
my $primer_fasta=$opt_p;
my $margin=(!defined($opt_m)?$MARGIN:$opt_m);
my $output_root=$opt_o;
my $estim_ampl_len=(!defined($opt_a)?0:$opt_a);

my $pattern_fasta_name="$output_root\.pattern";
my $search_output="$output_root\.search";
my $report_output="$output_root\.stats";
my $clear_range_output="$output_root\.clear";
my $no_primers_output="$output_root\.no_primers";

###############################################################################
# Generate pattern file based on primer.fasta file

print STDERR "Making Pattern File...\n";

open(PATTERN_FH, ">$pattern_fasta_name") || die "Could not open $pattern_fasta_name for writing.\n";
open(PRIMER_FH, "<$primer_fasta") || die "Could not open $primer_fasta for reading.\n";

my %primer_sequence_hash;

my ($defline, $prev_defline, $sequence);
while(<PRIMER_FH>){
        chomp;
        if(/^>/){
                $defline=$_;
                if($sequence ne ""){
                        process_primer_record($prev_defline, $sequence);
                        $sequence="";
                }
                $prev_defline=$defline;
        }else{
                $sequence.=$_;
        }
}
process_primer_record($prev_defline, $sequence);

close(PRIMER_FH);
close(PATTERN_FH);

#------------------------------------------------------------------------------

my %primer_params;

sub process_primer_record{
        my $defline = shift;
        my $sequence = shift;

	# Get id
	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Error parsing defline: $_\n";
	}

	# store primer sequence for future usage
	$primer_sequence_hash{$id}=$sequence;

	# Compute Error rate
	my $length=length($sequence);
	my $primer_error_rate;
	if(defined($opt_E)){
		$primer_error_rate=int($length*(1.0-$min_perc_id));
	}else{
		$primer_error_rate=$error_rate
	}

	# Compute actual percent identity
	my $actual_perc_id = ($length-$primer_error_rate)/$length;
	if($min_perc_id != $actual_perc_id){
		if(!defined($opt_E)){$min_perc_id="NA";}
		print STDERR "For $id: Target Identity: $min_perc_id, Applicable Identity: $actual_perc_id\n";
	}

	# Keep track of per primer parameters for later
	$primer_params{$id}{"EFF_PERC_ID"}=$actual_perc_id;
	$primer_params{$id}{"MAX_ERROR_RATE"}=$primer_error_rate;

	# Output FASTA
	print PATTERN_FH ">$id <mismatch=$primer_error_rate>\n";
	my $width=60;
        my $pos=0;
        do{
                my $out_width=($width>$length)?$length:$width;
                print PATTERN_FH substr($sequence, $pos, $width) . "\n";
                $pos+=$width;
                $length-=$width;
        }while($length>0);
}


###############################################################################
# Read in the read sequences that should be trimmed

my %read_length_hash;
my $total_input_sequences=0;

open(READS_FH, "<$reads_fasta") || die "Could not open $reads_fasta for reading.\n";

$sequence="";
while(<READS_FH>){
        chomp;
        if($_=~/^>/){
                $defline=$_;
                if($sequence ne ""){
                        process_reads_record($prev_defline, $sequence);
                        $sequence="";
                }
                $prev_defline=$defline;
        }else{
                $sequence.=$_;
        }
}
process_reads_record($prev_defline, $sequence);

close(READS_FH);

#------------------------------------------------------------------------------

sub process_reads_record{
	my $defline=shift;
	my $sequence=shift;

	$total_input_sequences++;
	if($defline=~/^>(\S+)/){
		my $id=$1;
		$read_length_hash{$id}=length($sequence);
	}else{
		die "Could not parse defline of read: $defline\n";
	}
}

###############################################################################
# Execute fuzznuc search, filter results, and produce output

# Setup fuzznuc parameters and execute
my $exec_string = join " ", (
		$FUZZ_NUC_BIN, 
		"-filter", 
		"-complement", "Yes",
		"-sequence", $reads_fasta,
		"-pattern", "@" . "$pattern_fasta_name"
		);
print STDERR "Executing: '$exec_string'\n";
print STDERR `$exec_string > $search_output`;

# Parse results from fuzznuc
my $results_hash=parse_fuzznuc($search_output);

# Filter results and generate filtering report
my ($filtered_hash, $report_hash_ref)=filter_results($results_hash, $margin, $estim_ampl_len);

# Output filter report
${$report_hash_ref}{"PARAM.ESTIM_AMPL_LEN"}=$estim_ampl_len;
${$report_hash_ref}{"PARAM.MARGIN"}=$margin;
foreach my $primer_name(keys %primer_params){
	${$report_hash_ref}{"PARAM.$primer_name.EFF_PERC_ID"}=$primer_params{$primer_name}{"EFF_PERC_ID"};
	${$report_hash_ref}{"PARAM.$primer_name.MAX_ERROR_RATE"}=$primer_params{$primer_name}{"MAX_ERROR_RATE"};
}
open(REPORT_FH, ">$report_output") || die "Could not open $report_output for writing.\n";
foreach my $key(sort keys %{$report_hash_ref}){
	print REPORT_FH "$key=${$report_hash_ref}{$key}\n";
}
close(REPORT_FH);

# Output clear coordinates
open(CLEARRANGE_FH, ">$clear_range_output") || die "Could not open $clear_range_output for writing.\n";
open(NOPRIMERS_FH, ">$no_primers_output") || die "Could not open $no_primers_output for writing.\n";
foreach my $seq_id (sort keys %read_length_hash){
	my ($begin, $end, $primer_info)=compute_clear_range(${$filtered_hash}{$seq_id}, $read_length_hash{$seq_id}, $margin, $estim_ampl_len);

	if($begin eq $end){
		print STDERR "$seq_id: Suspect primers locations detected.  Should be trashed.\n";
	}

	print CLEARRANGE_FH "$seq_id\t$begin\t$end\t$primer_info\n";

	if(($end-$begin)==$read_length_hash{$seq_id}){
		print NOPRIMERS_FH "$seq_id\n";
	}
}

close(CLEARRANGE_FH);
close(NOPRIMERS_FH);

###############################################################################


sub compute_clear_range{
	my $hit_list_ref=shift;
	my $read_len=shift;
	my $margin=shift;
	my $est_amp_len=shift;

	my ($clear_begin, $clear_end)=(0, $read_len);
	my ($primers_on_left, $primers_on_right)=(0,0);

	my @primer_info;

	my $extra_short_amp_padding=0;
	if(defined($est_amp_len) && $est_amp_len!=0){
		$extra_short_amp_padding=$read_len-$est_amp_len;
	}

	foreach my $hit_list_rec(@{$hit_list_ref}){
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$hit_list_rec};

		my $left=max($start, $end);
		my $right=min($start, $end);

		my $which_end;
	
		if(($left<$margin) && ($right>($read_len-($margin+$extra_short_amp_padding)))){
			# Trash the whole read
			$clear_begin=0;
			$clear_end=0;
			$which_end="amb";
		}elsif($left<$margin){
			$clear_begin=$left;
			$primers_on_left++;
			$which_end="five";
		}elsif($right>($read_len-($margin+$extra_short_amp_padding))){
			$clear_end=$right;
			$primers_on_right++;
			$which_end="three";
		}

		push @primer_info, "$primer_name/$mismatch_count/$which_end";
	}

	if($primers_on_left>1 || $primers_on_right>1){
		# Trash the whole read
		$clear_begin=0;
		$clear_end=0;
	}

	return($clear_begin, $clear_end, join ";", @primer_info);
}

###############################################################################

sub max{
	$a=shift;
	$b=shift;
	return(($a>$b)?$a:$b);
}
sub min{
	$a=shift;
	$b=shift;
	return(($a<$b)?$a:$b);
}

#------------------------------------------------------------------------------

sub filter_results{
	my $results_hash_ref=shift;
	my $margin=shift;
	my $est_amp_len=shift;
	my %info_hash;

	$info_hash{"MARGIN.5_PRIME_KEPT"}=0;
	$info_hash{"MARGIN.3_PRIME_KEPT"}=0;
	$info_hash{"MARGIN.FILTERED"}=0;
	$info_hash{"BESTHIT.KEPT"}=0;
	$info_hash{"BESTHIT.FILTERED"}=0;

	my $BINS=30;
	my $READ_LENGTH=1700;
	my $BIN_SIZE=$READ_LENGTH/$BINS;
	for(my $i=0; $i<$BINS; $i++){
		my $bin_str=sprintf("%04i", $i*$BIN_SIZE);
		$info_hash{"HISTOGRAM.VALUES_$bin_str"}=0;
	}

	my $total_sequences_w_hits=0;
	foreach my $seq_id(sort keys %{$results_hash_ref}){
		#print ">>>>>>>>$seq_id:\n";
		my $hits_arr_ref=${$results_hash_ref}{$seq_id};
		_filter_by_besthit($hits_arr_ref, \%info_hash);
		_filter_by_margin($hits_arr_ref, \%info_hash, $margin, $est_amp_len, $read_length_hash{$seq_id});
		_filter_inside_hits($hits_arr_ref, \%info_hash, $margin, $read_length_hash{$seq_id});
		_generate_histogram($hits_arr_ref, \%info_hash, $BIN_SIZE);
		${$results_hash_ref}{$seq_id}=$hits_arr_ref;
		$total_sequences_w_hits+=($#{$hits_arr_ref}>=0);
	}
	$info_hash{"TOTAL.NUM_SEQ_W_PRIMERS"}=$total_sequences_w_hits;
	$info_hash{"TOTAL.INPUT_SEQUENCES"}=$total_input_sequences;
	$info_hash{"TOTAL.WITHOUT_PRIMERS"}=$total_input_sequences-$total_sequences_w_hits;
	$info_hash{"TOTAL.PERC_W_PRIMERS"}=sprintf("%4.2f%%", 100.0*$total_sequences_w_hits/$total_input_sequences);

	# Put all kept hits into one array
	my @all_hits;
	foreach my $seq_id(sort keys %{$results_hash_ref}){
		foreach my $hit(@{${$results_hash_ref}{$seq_id}}){
			push @all_hits, $hit;
		}	
	}
	_profile_primer_mismatches(\@all_hits, \%info_hash, \%primer_sequence_hash);

	return($results_hash_ref, \%info_hash);
}

#-------------------------------------------------------------------------------

sub _generate_histogram{
	my $arr_ref=shift;
        my $info_hash_ref=shift;
	my $bin_size=shift;

        my @arr=@{$arr_ref};
	foreach my $rec_ref(@{$arr_ref}){
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
		my $midpoint=int(($start+$end)/2.0);
		my $bin=int($midpoint/$bin_size);
		my $bin_str=sprintf("%04i", $bin*$bin_size);
		${$info_hash_ref}{"HISTOGRAM.VALUES_$bin_str"}++;
	}
}

#------------------------------------------------------------------------------

sub _filter_by_besthit{
	my $arr_ref=shift;
	my $info_hash_ref=shift;

	my @arr=@{$arr_ref};
	my @kept_arr=();

	my $num_kept=0;
	my $num_filtered=0;

	#print "BEFORE:\n";
	#foreach my $rec_ref(@arr){
	#	my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
	#	print "\t$start, $end, $primer_name, $mismatch_count, $sequence\n";
	#}

	# Split primer hits by primer name
	my %primer_specific_hits;
	foreach my $hit(@arr){
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$hit};
		push @{$primer_specific_hits{$primer_name}}, $hit;
	}

	# Sort by decreasing mismatch score
	sub by_mismatch{
		my ($astart, $aend, $aprimer_name, $amismatch_count, $asequence)=@{$a};
		my ($bstart, $bend, $bprimer_name, $bmismatch_count, $bsequence)=@{$b};
		$amismatch_count <=> $bmismatch_count;
	}

	foreach my $primer_name(keys %primer_specific_hits){

		# Sort hits for each primer
		@{$primer_specific_hits{$primer_name}}=sort by_mismatch @{$primer_specific_hits{$primer_name}};

		# Keep record as good as best
		my $best_hit_mismatch_count=${${$primer_specific_hits{$primer_name}}[0]}[3];
		foreach my $rec_ref(@{$primer_specific_hits{$primer_name}}){
			my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
			if($mismatch_count==$best_hit_mismatch_count){
				push @kept_arr, $rec_ref;
				$num_kept++;
			}else{
				$num_filtered++;
			}
		}
	}
	
	#print "AFTER\n";
	#foreach my $rec_ref(@kept_arr){
	#	my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
	#	print "\t$start, $end, $primer_name, $mismatch_count, $sequence\n";
	#}

	${$info_hash_ref}{"BESTHIT.KEPT"}+=$num_kept;
	${$info_hash_ref}{"BESTHIT.FILTERED"}+=$num_filtered;

	@{$arr_ref}=@kept_arr;
}

#-------------------------------------------------------------------------------

sub _filter_by_margin{
	my $arr_ref=shift;
	my $info_hash_ref=shift;
	my $margin=shift;
	my $est_ampl_len=shift;
	my $read_length=shift;

	my @arr=@{$arr_ref};
	my @kept_arr=();

	my $num_kept_fiveprime=0;
	my $num_kept_threeprime=0;
	my $num_filtered=0;

	my $extra_short_amp_padding=0;
	if(defined($est_ampl_len) && $est_ampl_len!=0){
		$extra_short_amp_padding=$read_length-$est_ampl_len;
	}

	foreach my $rec_ref(@arr){
		my @hit_info=@{$rec_ref};
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@hit_info;
		if(max($start, $end)<$margin){
			push @kept_arr, $rec_ref;
			$num_kept_fiveprime++
		}elsif(min($start, $end)>=($read_length-($margin+$extra_short_amp_padding))){
			push @kept_arr, $rec_ref;
			$num_kept_threeprime++;
		}else{
			$num_filtered++;
		}
	}

	${$info_hash_ref}{"MARGIN.5_PRIME_KEPT"}+=$num_kept_fiveprime;
	${$info_hash_ref}{"MARGIN.3_PRIME_KEPT"}+=$num_kept_threeprime;
	${$info_hash_ref}{"MARGIN.FILTERED"}+=$num_filtered;

	@{$arr_ref}=@kept_arr;
}

#-------------------------------------------------------------------------------

sub _filter_inside_hits{
	my $arr_ref=shift;
	my $info_hash_ref=shift;
	my $margin=shift;
	my $read_length=shift;

	my @arr=@{$arr_ref};
	my @kept_arr=();

	my $num_filtered=0;
	my $num_kept=0;

	my %primer_name_hash;
	my %primer_info_hash;
	my %primer_edge_dist_hash;

	#print "readlen=$read_length\n";
	#foreach my $rec_ref(@arr){
	#	my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
	#	print "BEFORE:\t$start, $end, $primer_name, $mismatch_count, $sequence\n";
	#}
	
	#print "------------------------\n";

	foreach my $rec_ref(@arr){
		my @hit_info=@{$rec_ref};
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@hit_info;

		# Split the hits by primer name
		push @{$primer_name_hash{$primer_name}}, $rec_ref;
		# Key the hits by position
		$primer_info_hash{"$start#$end"}=$rec_ref;

		# Compute distance from edge
		my $from_5prime=max($start,$end);
		my $from_3prime=$read_length-min($start,$end);
		my $closest_to_edge=min($from_5prime, $from_3prime);

		# Keep track of edge distance to record association
		push @{$primer_edge_dist_hash{$primer_name}}, "$closest_to_edge,$start#$end";
	}

	sub by_edge_dist{
		my ($adist, $arecord)=split /,/, $a;
		my ($bdist, $brecord)=split /,/, $b;
		$adist<=>$bdist;
	}

	foreach my $primer_name(keys %primer_edge_dist_hash){
		# Sort the records for this primer by distance from edge
		my @sorted_by_edge_dist=sort by_edge_dist @{$primer_edge_dist_hash{$primer_name}};
		my ($dist, $rec_key)=split /,/, $sorted_by_edge_dist[0];
		push @kept_arr, $primer_info_hash{$rec_key};
		$num_kept++;
		$num_filtered+=$#sorted_by_edge_dist;
	}

	#foreach my $rec_ref(@kept_arr){
	#	my ($start, $end, $primer_name, $mismatch_count, $sequence)=@{$rec_ref};
	#	print "AFTER:\t$start, $end, $primer_name, $mismatch_count, $sequence\n";
	#}

	${$info_hash_ref}{"INSIDE.FILTERED"}+=$num_filtered;
	${$info_hash_ref}{"INSIDE.KEPT"}+=$num_kept;

	@{$arr_ref}=@kept_arr;
}

#------------------------------------------------------------------------------

sub _profile_primer_mismatches{
	my $arr_ref=shift;
	my $info_hash_ref=shift;
	my $primer_sequence_hash_ref=shift;
	my $primer_sequence_hash=%{$primer_sequence_hash_ref};

	my %recs_by_primer;

	# Split hits by primer name, keep only the sequence.
	foreach my $rec_ref(@{$arr_ref}){
		my @hit_info=@{$rec_ref};
		my ($start, $end, $primer_name, $mismatch_count, $sequence)=@hit_info;
		push @{$recs_by_primer{$primer_name}}, $sequence;
	}

	# Count up positions/errors
	my %primer_profile;
	my %total_reads_for_primer;
	foreach my $primer_name(keys %recs_by_primer){
		foreach my $match_sequence(@{$recs_by_primer{$primer_name}}){
			$match_sequence=uc($match_sequence);
			my @match_seq_arr=split //, $match_sequence;
			for(my $i=0; $i<=$#match_seq_arr; $i++){
				$primer_profile{$primer_name}[$i]{$match_seq_arr[$i]}++;
			}
			$total_reads_for_primer{$primer_name}++;
		}
	}

	# Save ouptut
	foreach my $primer_name(keys %recs_by_primer){

		${$info_hash_ref}{"PROFILE\.$primer_name\.TOTALS"}=$total_reads_for_primer{$primer_name};

		my @primer_seq_arr=split //, $primer_sequence_hash{$primer_name};
		for(my $i=0; $i<=$#primer_seq_arr; $i++){
			my @prof_arr;
			foreach my $nuc("A", "T", "G", "C", "N", "X"){
				my $ind_prof=sprintf("%3.1f", 100.0*$primer_profile{$primer_name}[$i]{$nuc}/$total_reads_for_primer{$primer_name});
				if(!defined($ind_prof)){
					$ind_prof=0;
				}		
				push @prof_arr, "$nuc=$ind_prof";
			}
			my $prof_string=join " ", @prof_arr;
			my $nuc_offset=sprintf("%02i", $i);
			${$info_hash_ref}{"PROFILE\.$primer_name\.$nuc_offset\.$primer_seq_arr[$i]"}.="{$prof_string}";
		}
	}


}

#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------

sub parse_fuzznuc{

	my $search_output=shift;
	my %results_hash;

	open(FZNC_FH, "<$search_output") || die "Could not open $search_output\n";

	my $current_seq_id;
	while(<FZNC_FH>){
		chomp;
		if($_=~/^# Sequence: (\S+) /){
			$current_seq_id=$1;
		}elsif($_=~/^#/ || $_=~/^$/){
			# Noop
		}elsif($_=~/^  Start\s+End\s+Pattern_name\s+Mismatch\s+Sequence/){
			# Noop
		}else{
			$_=~s/^\s+//;
			my ($start, $end, $primer_name, $mismatch_count, $sequence)=split /\s+/, $_;
			$start--; #convert to 0-space based coordinates 	
			$mismatch_count=($mismatch_count eq ".")?0:int($mismatch_count);
			#print ((join "\t", ($start, $end, $primer_name, $mismatch_count, $sequence)) . "\n");
			my @record= ($start, $end, $primer_name, $mismatch_count, $sequence);
			push @{$results_hash{$current_seq_id}}, \@record;
		}

	}

	close(FZNC_FH);
	return(\%results_hash);
}
