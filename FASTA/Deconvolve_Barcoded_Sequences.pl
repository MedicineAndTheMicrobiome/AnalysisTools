#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use Sys::Hostname;
use File::Basename;
use Cwd;
use vars qw($opt_f $opt_b $opt_o $opt_p $opt_c $opt_r $opt_d $opt_a $opt_s);
getopts("f:b:o:p:c:r:d:as");

# Need at least EMBOSS 5.0, 4 is broken
my $FUZZ_NUC_BIN="/usr/local/packages/EMBOSS-5.0.0/bin/fuzznuc";
my $fuzznuc_found_str=`which $FUZZ_NUC_BIN`;
if($fuzznuc_found_str=~/Command not found/){
	print STDERR "$FUZZ_NUC_BIN not found.\n";
}

# Defaults
my $GOODHITS_FNAME="good_hits";
my $NOHITS_FNAME="no_hits";
my $ERRORHITS_FNAME="error_hits";
my $BADHITS_FNAME="bad_hits";
my $STATISTICS_FNAME="stats";

my $SPLIT_DIR="deconvolved";

my $BARCODE_ERROR_RATE=1;
my $PRIMER_ERROR_RATE=6;
my $ADAPTER_ERROR_RATE=0;

my $MAX_BARCODE_PRIMER_DISTANCE=3;
my $MAX_BARCODE_END_DISTANCE=2;
my $MAX_ADAPTER_END_DISTANCE=1;

my $UNBARCODE_FILENAME="Unknown";
my $FASTA_EXTENSION="fasta";
my $PATPREFIX="OligoPat";

# Adapter A
# PrimerA1 : CCATCTCATCCCTGCGTGTCTCCGACTCAG
# PrimerA'1:                   TCTCCGACTCAG

# Adapter B
# PrimerB1 : CCTATCCCCTGTGTGCCTTGGCAGTCTCAG
# PrimerB'1:                   TGGCAGTCTCAG

my $ADAPTER_A_DS_SEQ="TCTCCGACTCAG";
my $ADAPTER_B_DS_SEQ="TGGCAGTCTCAG";

my $usage="
	
	$0 

   Required Parameters:
	-f <reads fasta sequence file>
		Contains barcoded reads you want to deconvolve (FASTA).

	-b <barcode to sample id map file>
		Contains list of barcode sequence to sample id mapping (one sequence per line, not FASTA).
		Only the first column is utilized.

	-o <output directory name>
		Output directory where deconvolved data will sit.  Directory will be created if it does not exist.

   Optional Parameters:
	Recommended:
	    [-p <primer sequence list, DEFAULT, none used.> ]
		    Contains list of primer sequences that the barcodes are connected to (one sequence per line, not FASTA).

	Tolerance-related:
	    [-c <barcode error rate, DEFAULT=$BARCODE_ERROR_RATE >]
		    Number of bases difference to allow in barcode.
	    [-r <primer error rate, DEFAULT=$PRIMER_ERROR_RATE >]
		    Number of bases difference to allow in primer.
	    [-d <barcode-to-primer max distance, DEFAULT=$MAX_BARCODE_PRIMER_DISTANCE >]

	Other:
	    [-a <flag: Do not search and trim out double stranded segment of A and B adapter: $ADAPTER_A_DS_SEQ/$ADAPTER_B_DS_SEQ>]
	    [-s <flag: Save intermediate search and pattern files>]

Example usage:

    Deconvolve_Barcoded_Sequences.pl -f FUXKOPD02.fasta -b barcode_sample.map -o FUXKOPD02_deconv -p primers


This program will search through the input reads and try to sort out which reads should be associated with what barcode.  The output
will consist of several text files:

$GOODHITS_FNAME : Contains reads that are in the correct orientation and position
$NOHITS_FNAME : Contains reads that have no barcoded detected in them whatsoever
$ERRORHITS_FNAME : Contains reads that have barcodes but for some additional reason are problematic, for example multi-barcoded reads
$BADHITS_FNAME : Contains reads where the barcode was found, but not in the correct orientation or position
$STATISTICS_FNAME : Contains the percentages of the above classifications


Example output for $GOODHITS_FNAME, $ERRORHITS_FNAME, and $BADHITS_FNAME:

    read_id         barcode       start  end  ori  mismatch  clr_begin clr_end  read_len  barcode_end  primer_hit      adapter_trim
    FUXKOPD01AH8IF  CGATACTAGACG  0      12   1    0         12        126      126       5';0         12;32;1;1;0     0-126 (0,0)
    FUXKOPD01AHO6W  CGCTGAGTACTC  0      12   1    0         12        215      215       5';0         12;29;1;0;0     0-215 (0,0)
    FUXKOPD01AI23J  CGTACGTCTGTG  96     108  -1   0         0         96       108       3';0         79;96;-1;2;0    0-108 (0,0)
    FUXKOPD01AI9EH  CGAGTATCGATG  0      12   1    0         12        79       79        5';0         12;29;1;2;0     0-79 (0,0)

	The barcode_end should be interpreted:
		<end>;<distance>
		    <end>: Whether barcode was detected on 5' or 3' end of read
		    <distance>:  How far away from the end, the barcoded was detected

	The primer_hit should be interpreted:
		<begin>;<end>;<orientation>;<primer_id>;<mismatches>

		<primer_id> is the ordinal number of the primer once it has been sorted alphabetically start from 0

	The adapter_trim should be interpreted:
		<begin>-<end> (<num_begin_trims>,<num_end_trims>)
	
		<begin>,<end> represent the range in the read where the adapter was detected
		<num_begin_trims>,<num_end_trims> represents the number of times the adapter was detected.

	The clr_begin and clr_end can be used for downstream sequence trimming/masking

";


if(
	!defined($opt_f) ||
	!defined($opt_b) ||
	!defined($opt_o)
){
	die $usage;
}

my $input_fasta_file=$opt_f;
my $barcode_map_file=$opt_b;
my $output_dir=$opt_o;
my $primers_map_file=$opt_p;
my $trim_a_adapter=(!defined($opt_a));
my $save_tmp_files=defined($opt_s);

print STDERR "Save temp files?: $save_tmp_files\n";

my $barcode_error_rate=$BARCODE_ERROR_RATE;
my $primer_error_rate=$PRIMER_ERROR_RATE;
my $adapter_error_rate=$ADAPTER_ERROR_RATE;
my $barcode_primer_max_distance=$MAX_BARCODE_PRIMER_DISTANCE;

if(defined($opt_c)){
	$barcode_error_rate=$opt_c;
}

if(defined($opt_r)){
	$primer_error_rate=$opt_r;
}

if(defined($opt_d)){
	$barcode_primer_max_distance=$opt_d;
}

if(substr($output_dir,0,1) ne "/"){
	$output_dir= cwd() . "/$output_dir";
}

my $hostname=hostname();
my $pid=$$;
my ($progname)=fileparse($0);
my $temp_filename_root="$output_dir/$progname\.$hostname\.$pid";

my $BARCODES_PAT_EXT="barcodes.pat";
my $PRIMERS_PAT_EXT="primers.pat";
my $ADAPTERS_PAT_EXT="adapters.pat";

my $BARCODES_SEARCH_RESULT_EXT="barcodes.search.out";
my $PRIMERS_SEARCH_RESULT_EXT="primers.search.out";
my $ADAPTERS_SEARCH_RESULT_EXT="adapters.search.out";

###############################################################################
# Generate pattern file based on primer.fasta file

if(!(-e $output_dir)){
	print STDERR "Creating $output_dir\n";
	`mkdir $output_dir`;
}else{
	print STDERR "Using existing $output_dir\n";
}

# Load barcodes
print STDERR "Loading Barcodes...\n";
my $barcodes_hash_ref=load_oligos($barcode_map_file);
my @barcodes_arr= sort keys %{$barcodes_hash_ref};
my $num_barcodes=$#barcodes_arr+1;
print STDERR "Num barcodes loaded: $num_barcodes\n";

# Load primers
my $primers_hash_ref;
my @primers_arr;
my $num_primers;
if(defined($primers_map_file)){
	print STDERR "Primer file specified, so loading these now...\n";
	$primers_hash_ref=load_oligos($primers_map_file);
	@primers_arr= sort keys %{$primers_hash_ref};
	$num_primers=$#primers_arr + 1;
	print STDERR "Num primers loaded: $num_primers\n";
}

# Load reads
my $read_id_hash_ref=load_reads($input_fasta_file);
my @reads_arr=keys %{$read_id_hash_ref};
my $num_reads=$#reads_arr+1;
print STDERR "Num Reads to deconvolve: $num_reads\n";

# Setup array of adapter sequences
my @adapter_arr=($ADAPTER_A_DS_SEQ, $ADAPTER_B_DS_SEQ);

# Setup patterns file for barcodes and primers
make_pattern_file(\@adapter_arr, $adapter_error_rate, "$temp_filename_root\.$ADAPTERS_PAT_EXT");
make_pattern_file(\@barcodes_arr, $barcode_error_rate, "$temp_filename_root\.$BARCODES_PAT_EXT");
if(defined($primers_map_file)){
	make_pattern_file(\@primers_arr, $primer_error_rate, "$temp_filename_root\.$PRIMERS_PAT_EXT");
}

# Run fuzznuc
if(defined($trim_a_adapter)){
	run_fuzznuc("$temp_filename_root\.$ADAPTERS_PAT_EXT", $input_fasta_file, "$temp_filename_root\.$ADAPTERS_SEARCH_RESULT_EXT");
}
run_fuzznuc("$temp_filename_root\.$BARCODES_PAT_EXT", $input_fasta_file, "$temp_filename_root\.$BARCODES_SEARCH_RESULT_EXT");
if(defined($primers_map_file)){
	run_fuzznuc("$temp_filename_root\.$PRIMERS_PAT_EXT" , $input_fasta_file, "$temp_filename_root\.$PRIMERS_SEARCH_RESULT_EXT");
}

# Parse fuzznuc results
my $adapters_hits_results_hash_ref=undef;
if(defined($trim_a_adapter)){
	my $no_keep;
	($adapters_hits_results_hash_ref, $no_keep)=parse_fuzznuc("$temp_filename_root\.$ADAPTERS_SEARCH_RESULT_EXT");
}
my ($barcode_hits_results_hash_ref, $read_len_hash_ref)=parse_fuzznuc("$temp_filename_root\.$BARCODES_SEARCH_RESULT_EXT");
my $primers_hits_results_hash_ref=undef;
if(defined($primers_map_file)){
	my $no_keep;
	($primers_hits_results_hash_ref, $no_keep)=parse_fuzznuc("$temp_filename_root\.$PRIMERS_SEARCH_RESULT_EXT");
}

#foreach my $read_id(keys %{$barcode_hits_results_hash_ref}){
#	print "$read_id\n";
#	foreach my $align (@{${$barcode_hits_results_hash_ref}{$read_id}}){
#		my $align_str=join " ", @{$align};
#		print "\t$align_str\n";
#	}
#}

# Filter barcode hits
my $primer_required=0;
my ($good_hits_ref, $bad_hits_ref, $multi_hits_ref, $num_reads_w_primersupport, $num_reads_w_endsupport, $hit_info_hash_ref)=
	filter_hits(	$barcode_hits_results_hash_ref, $primers_hits_results_hash_ref, $adapters_hits_results_hash_ref,
			\@barcodes_arr, \@primers_arr, \@adapter_arr,
			$read_len_hash_ref, $primer_required);

# Recover double barcoded reads with the same barcode
recover_double_barcoded_hits($multi_hits_ref, $good_hits_ref, $hit_info_hash_ref);

# Output 
my $num_nohits=dump_nohits(\@reads_arr, $barcode_hits_results_hash_ref, "$output_dir/$NOHITS_FNAME");

# Output Results
dump_results($good_hits_ref, \@barcodes_arr, $read_len_hash_ref, $hit_info_hash_ref, "$output_dir/$GOODHITS_FNAME");
dump_results($bad_hits_ref, \@barcodes_arr, $read_len_hash_ref, $hit_info_hash_ref, "$output_dir/$BADHITS_FNAME");
dump_results($multi_hits_ref, \@barcodes_arr, $read_len_hash_ref, $hit_info_hash_ref, "$output_dir/$ERRORHITS_FNAME");

# Output Stats
open(STAT_FH, ">$output_dir/$STATISTICS_FNAME") || die "Could not open $output_dir/$STATISTICS_FNAME";

print STAT_FH "\n\nNum Reads: $num_reads\n";
my $perc;

$perc=sprintf("%3.2f", 100*($num_reads_w_primersupport/$num_reads));
print STAT_FH "Num Reads with Primer Support: $num_reads_w_primersupport ($perc%)\n";
$perc=sprintf("%3.2f", 100*($num_reads_w_endsupport/$num_reads));
print STAT_FH "Num Reads with End Support: $num_reads_w_endsupport ($perc%)\n\n";

my @arr;
my $count;

@arr=keys %{$good_hits_ref};
$count=$#arr+1;
$perc=sprintf("%3.2f", 100*$count/$num_reads);
print STAT_FH "Num Good Reads: " . $count . " ($perc%)\n";

@arr=keys %{$multi_hits_ref};
$count=$#arr+1;
$perc=sprintf("%3.2f", 100*$count/$num_reads);
print STAT_FH "Num Multibarcoded Reads: " . $count . " ($perc%)\n";

@arr=keys %{$bad_hits_ref};
$count=$#arr+1;
$perc=sprintf("%3.2f", 100*$count/$num_reads);
print STAT_FH "Num Bad Hit Reads: " . $count . " ($perc%)\n";

$count=$num_nohits;
$perc=sprintf("%3.2f", 100*$count/$num_reads);
print STAT_FH "Num No Hit Reads: " . $count . " ($perc%)\n";

close(STAT_FH);

if(!$save_tmp_files){
	`rm $temp_filename_root\.$BARCODES_PAT_EXT`;
	`rm $temp_filename_root\.$BARCODES_SEARCH_RESULT_EXT`;

	if($trim_a_adapter){
		`rm $temp_filename_root\.$ADAPTERS_PAT_EXT`;
		`rm $temp_filename_root\.$ADAPTERS_SEARCH_RESULT_EXT`;
	}

	if(defined($primers_map_file)){
		`rm $temp_filename_root\.$PRIMERS_PAT_EXT`;
		`rm $temp_filename_root\.$PRIMERS_SEARCH_RESULT_EXT`;
	}
}

###############################################################################

sub recover_double_barcoded_hits{
	my $multi_hits_ref=shift;
	my $good_hits_ref=shift;
	my $hit_info_hash_ref=shift;
	my $num_recovered=0;

	print STDERR "Recovering multi-barcodes...\n";

	foreach my $id(keys %{$multi_hits_ref}){

		#print STDERR "$id\n";
		my $hits_ref=${$multi_hits_ref}{$id};

		if($#{$hits_ref}==1){ 		# If there are exactly two barcode hits

			my $hit_a=${$hits_ref}[0];
			my $hit_b=${$hits_ref}[1];

			my ($astart, $aend, $aori, $aoligo_name, $amismatch_count)=split /;/, $hit_a;
			my ($bstart, $bend, $bori, $boligo_name, $bmismatch_count)=split /;/, $hit_b;

			if(($aoligo_name eq $boligo_name) &&
			   ($aori ne $bori)){
				#print STDERR "Moving read from multi to good hit.\n";
				
				my $hit_comb=join "/", ($hit_a, $hit_b);

				# Insert combined record into good hits list
				push @{${$good_hits_ref}{$id}}, $hit_comb;
	
				# Remove records from multi hits list
				delete ${$multi_hits_ref}{$id};

				# Count recoveries
				$num_recovered++;
			}

		}
				
	}
	print STDERR "Num double barcoded hits recovered: $num_recovered\n";

}

#------------------------------------------------------------------------------

sub dump_nohits{
	my $reads_arr_ref=shift;
	my $barcode_hits_results_hash_ref=shift;
	my $outputfname=shift;
	my $no_hits=0;

	open(FH, ">$outputfname") || die "Could not open $outputfname for writing no hits file.\n";

	foreach my $read_id(@{$reads_arr_ref}){
		if(!defined(${$barcode_hits_results_hash_ref}{$read_id})){
			print FH "$read_id\n";
			$no_hits++;
		}
	}

	close(FH);

	return($no_hits);
}

#------------------------------------------------------------------------------

sub dump_results{
	my $hits_ref=shift;
	my $barcodes_arr_ref=shift;
	my $read_length_hash_ref=shift;
	my $hit_info_hash_ref=shift;
	my $fname=shift;

	open(FH, ">$fname") || die "Could not open $fname for writing.\n";

	my $title_string= join "\t",
		("read_id", "barcode", "start", "end", "ori", "mismatch", 
		"clr_begin", "clr_end", "read_len", "barcode_end", "primer_hit", "adapter_trim");
	print FH "$title_string\n";
	
	foreach my $id(sort keys %{$hits_ref}){
		
		# Get the read length
		my $read_len=${$read_length_hash_ref}{$id};

		# Pull out the evidence information
		my $evidence_ref=${$hit_info_hash_ref}{$id};
		my $barcode_evidence=${$evidence_ref}{"BARCODE"};
		my $primer_evidence=${$evidence_ref}{"PRIMER"};
		my $adapter_evidence=${$evidence_ref}{"ADAPTER"};

		my $adapter_info=$adapter_evidence;	# This is based on the read, not on a barcode hit

		if($#{${$hits_ref}{$id}}==-1){
			print FH "$id\n";
		}else{

			foreach my $barcode_hit_rec (@{${$hits_ref}{$id}}){

				my @barcode_hits_arr=split /\//, $barcode_hit_rec;

				my $barcode;
				my @start_arr;
				my @end_arr;
				my @ori_arr;
				my @mismatch_count_arr;
				my @barcode_info_arr;
				my @primer_info_arr;

				@barcode_hits_arr=sort hits_by_begin_ascending @barcode_hits_arr;

				foreach my $barcode_hit(@barcode_hits_arr){

					my ($start, $end, $ori, $oligo_name, $mismatch_count)=split /;/, $barcode_hit;

					# Get the barcode sequence from the oligo_name
					$barcode=${$barcodes_arr_ref}[$oligo_name];

					push @start_arr, $start;
					push @end_arr, $end;
					push @ori_arr, $ori;
					push @mismatch_count_arr, $mismatch_count;

					# Get info on barcode hit
					push @barcode_info_arr, ${$barcode_evidence}{$barcode_hit};	# eg. 5';0

					# Get primer info on barcode hit
					my $primer_info;
					if(defined(${$primer_evidence}{$barcode_hit})){
						$primer_info=${$primer_evidence}{$barcode_hit};
					}else{
						$primer_info="no_primer_hits";
					}
					push @primer_info_arr, $primer_info;

				}


				# Compute the clear range
				my ($clr_begin, $clr_end);
				if($#barcode_hits_arr==0){
					if($ori_arr[0]==1){
						$clr_begin=$end_arr[0];
						$clr_end=$read_len;
					}elsif($ori_arr[0]==-1){
						$clr_begin=0;
						$clr_end=$start_arr[0];
					}else{
						die "Unknown orientation: $ori_arr[0]\n";
					}
				}else{
					$clr_begin=$end_arr[0];
					$clr_end=$start_arr[1];
				}

				# Output the information
				my $out_str=join "\t", (
					$id, $barcode, 
					(join "/", @start_arr),
					(join "/", @end_arr),
					(join "/", @ori_arr),
					(join "/", @mismatch_count_arr),
					$clr_begin, $clr_end, $read_len,
					(join "/", @barcode_info_arr), 
					(join "/", @primer_info_arr), 
					$adapter_info
				);
				print FH "$out_str\n";
			}
		}
	}


	close(FH);
}

###############################################################################

sub filter_hits{
	my $barcode_hits_results_hash_ref=shift;
	my $primer_hits_results_hash_ref=shift;
	my $adapters_hits_results_hash_ref=shift;

	my $barcodes_arr_ref=shift;
	my $primers_arr_ref=shift;
	my $adapters_arr_ref=shift;

	my $read_len_hash_ref=shift;
	my $primer_required=shift;
	
	my %good_hits_hash;
	my %no_hits_hash;
	my %error_hits_hash;
	my %hit_info_hash;
	
	my @barcode_lengths;
	my @primer_lengths;

	# Compute barcode lengths
	for(my $i=0; $i<=$#{$barcodes_arr_ref}; $i++){
		$barcode_lengths[$i]=length(${$barcodes_arr_ref}[$i]);
	}
	# Compute primer lengths
	if(defined($primers_arr_ref)){
		for(my $i=0; $i<=$#{$primers_arr_ref}; $i++){
			$primer_lengths[$i]=length(${$primers_arr_ref}[$i]);
		}
	}

	my $num_reads=0;
	my $num_reads_w_primersupport=0;
	my $num_reads_w_endsupport=0;

	foreach my $id (keys %{$barcode_hits_results_hash_ref}){

		#print STDERR "read_id: $id\n";
		$num_reads++;
		
		my $barcode_hits_ref=${$barcode_hits_results_hash_ref}{$id};
		my $primer_hits_ref =${$primer_hits_results_hash_ref}{$id};
		my $adapter_hits_ref=${$adapters_hits_results_hash_ref}{$id};

		# pick the barcode that is upstream to the primer
		my ($primer_supported_barcodes_ref, $primer_support);
		($primer_supported_barcodes_ref, $primer_support, ${$hit_info_hash{$id}}{"PRIMER"})=
			filter_by_primer_proximity($barcode_hits_ref, $primer_hits_ref);
		$num_reads_w_primersupport+=($primer_support>0);

		# trim off the adapters
		my ($clr_begin, $clr_end);
		($clr_begin, $clr_end, ${$hit_info_hash{$id}}{"ADAPTER"})=
			trim_adapters($adapter_hits_ref, $MAX_ADAPTER_END_DISTANCE, ${$read_len_hash_ref}{$id});
		
		# pick the barcode that is up against the 5' or 3' end and in the right orientation
		my ($end_supported_barcode_ref, $end_support);
		($end_supported_barcode_ref, $end_support, ${$hit_info_hash{$id}}{"BARCODE"})=
			filter_by_end_proximity($barcode_hits_ref, $MAX_BARCODE_END_DISTANCE, $clr_begin, $clr_end);
		$num_reads_w_endsupport+=($end_support>0);

		# If there are no good barcodes or more than one, then track it.
		if($end_support==1 && $primer_support<=1){
			$good_hits_hash{$id}=$end_supported_barcode_ref;	
		}elsif($end_support==0){
			$no_hits_hash{$id}=$barcode_hits_ref;
		}else{
			my $ORed_hits_ref=combine_hits_OR($primer_supported_barcodes_ref, $end_supported_barcode_ref);
			$error_hits_hash{$id}=$ORed_hits_ref;
		}

	}

	return(\%good_hits_hash, \%no_hits_hash, \%error_hits_hash, $num_reads_w_primersupport, $num_reads_w_endsupport, \%hit_info_hash);
}


###############################################################################

sub combine_hits_OR{
	my $hits_ref1=shift;
	my $hits_ref2=shift;
	my @combined_arr;

	my %hit_hash;

	foreach my $hit(@{$hits_ref1}){
		$hit_hash{$hit}=1;
	}

	foreach my $hit(@{$hits_ref2}){
		$hit_hash{$hit}=1;
	}

	@combined_arr=keys %hit_hash;

	return(\@combined_arr);

}

###############################################################################

sub filter_by_primer_proximity{
	my $barcode_hits_ref=shift;
	my $primer_hits_ref=shift;
	my $primer_required=shift;
	my @kept_hits;

	my $primer_support=0;
	my %unique_barcode_hash;
	my %hit_info_hash;		# stores barcode_hit to primer_hit relationship
	
	foreach my $barcode_hit(@{$barcode_hits_ref}){

		my ($bstart, $bend, $bori, $boligo_name, $bmismatch_count)=split /;/, $barcode_hit;
		#print STDERR "\tbarcode=($bstart, $bend, $bori, $boligo_name, $bmismatch_count) \n";

		if(!$primer_required && $#{$primer_hits_ref}==-1){
			push @kept_hits, $barcode_hit;
		}else{
			foreach my $primer_hit(@{$primer_hits_ref}){

				my ($pstart, $pend, $pori, $poligo_name, $pmismatch_count)=split /;/, $primer_hit;

				#print STDERR "\t\tprimer=($pstart, $pend, $pori, $poligo_name, $pmismatch_count)";

				if($bori == $pori){

					my $primer_to_barcode_distance;
					if($bori==1){
						$primer_to_barcode_distance=($pstart-$bend);
					}elsif($bori==-1){
						$primer_to_barcode_distance=($bstart-$pend);
					}else{
						die "unknown orientation: $bori\n";
					}

					if(abs($primer_to_barcode_distance)<$MAX_BARCODE_PRIMER_DISTANCE){	

						# Make sure a barcode is only supported once per primer hit.
						if(!defined($unique_barcode_hash{$barcode_hit})){
							push @kept_hits, $barcode_hit;
							$unique_barcode_hash{$barcode_hit}=1;
							$primer_support++;
							#print STDERR " Great!";
							$hit_info_hash{$barcode_hit}=$primer_hit;
						}
					}
				}
				#print STDERR "\n";
			}
		}
	}

	return(\@kept_hits, $primer_support, \%hit_info_hash);
}

#------------------------------------------------------------------------------

sub hits_by_begin_ascending{
	my ($astart, $aend, $aori, $aoligo_name, $amismatch_count)=split /;/, $a;
	my ($bstart, $bend, $bori, $boligo_nbme, $bmismbtch_count)=split /;/, $b;
	$astart<=>$bstart;
}

sub hits_by_end_descending{
	my ($astart, $aend, $aori, $aoligo_name, $amismatch_count)=split /;/, $a;
	my ($bstart, $bend, $bori, $boligo_nbme, $bmismbtch_count)=split /;/, $b;
	$aend<=>$bend;
}

sub split_and_sort_hits_by_orientation{
	my $hits_ref=shift;
	my @forward_hits;
	my @reverse_hits;

	# Split by orientation
	foreach my $hit(@{$hits_ref}){
		my ($start, $end, $ori, $oligo_name, $mismatch_count)=split /;/, $hit;
		if($ori==1){
			push @forward_hits, $hit;
		}elsif($ori==-1){
			push @reverse_hits, $hit;
		}else{
			die "Unknown orientation: $ori\n";
		}
	}

	# Sort forwards by begin
	@forward_hits=sort hits_by_begin_ascending @forward_hits;
	
	# Sort reverse by end in reverse
	@reverse_hits=sort hits_by_end_descending @reverse_hits;

	return(\@forward_hits, \@reverse_hits);

}

#------------------------------------------------------------------------------

sub trim_adapters{
	my $adapter_hits_ref=shift;
	my $max_adapter_end_distance=shift;
	my $read_length=shift;
	my $trim_info;				# Stores clear range and number of trims made

	# "Trim" off the adapter as many times as we need to.
	my ($five_prime_adapter_count, $three_prime_adapter_count)=(0, 0);
	my ($five_prime_begin, $three_prime_end)=(0, $read_length);

	my ($for_adapt_hits_ref, $rev_adapt_hits_ref)=split_and_sort_hits_by_orientation($adapter_hits_ref);

	foreach my $adapter_hit(@{$for_adapt_hits_ref}){
		my ($start, $end, $ori, $oligo_name, $mismatch_count)=split /;/, $adapter_hit;
		if(($start-$five_prime_begin)<=$max_adapter_end_distance){
			$five_prime_begin=$end;
			$five_prime_adapter_count++;
		}
	}

	foreach my $adapter_hit(@{$rev_adapt_hits_ref}){
		my ($start, $end, $ori, $oligo_name, $mismatch_count)=split /;/, $adapter_hit;
		if(($three_prime_end-$end)<=$max_adapter_end_distance){
			$three_prime_end=$start;
			$three_prime_adapter_count++;
		}
	}

	$trim_info="$five_prime_begin-$three_prime_end ($five_prime_adapter_count,$three_prime_adapter_count)";

	return($five_prime_begin, $three_prime_end, $trim_info);
}

#------------------------------------------------------------------------------

sub filter_by_end_proximity{
	my $barcode_hits_ref=shift;
	my $max_barcode_end_distance=shift;
	my $clr_begin=shift;
	my $clr_end=shift;

	my @kept_hits;
	my $end_support=0;
	my %hit_info_hash;

	# Determine if barcode is at end, based on clear coordinates
	foreach my $barcode_hit(@{$barcode_hits_ref}){

		my ($bstart, $bend, $bori, $boligo_name, $bmismatch_count)=split /;/, $barcode_hit;

		my $dist;
		if($bori==1){
			$dist=$bstart-$clr_begin;
			if($dist<$max_barcode_end_distance){
				push @kept_hits, $barcode_hit;
				$end_support++;
				$hit_info_hash{$barcode_hit}="5';$dist";
			}
		}elsif($bori==-1){
			$dist=$clr_end-$bend;
			if($dist<$max_barcode_end_distance){
				push @kept_hits, $barcode_hit;
				$end_support++;
				$hit_info_hash{$barcode_hit}="3';$dist";
			}
		}else{
			die "Unknown orientation: $bori\n";
		}
	}

	return(\@kept_hits, $end_support, \%hit_info_hash);
}


###############################################################################

sub run_fuzznuc{
	my $pattern_fasta_name=shift;
	my $reads_fasta=shift;
	my $search_output=shift;

	# Setup fuzznuc parameters and execute
	my $exec_string = join " ", (
			$FUZZ_NUC_BIN, 
			"-filter", 
			"-complement", "Yes",
			"-sequence", $reads_fasta,
			"-pattern", "@" . "$pattern_fasta_name"
			);
	print STDERR "Executing: '$exec_string'\n";

	# filter out the # $PARPREFIX because useful and takes up a lot of disk space in the output
	print STDERR `$exec_string | grep -v \"# $PATPREFIX\"  >  $search_output`;
}


###############################################################################

sub make_pattern_file{
	my $oligo_arr_ref=shift;
	my $error_rate=shift;
	my $output_file=shift;

	open(PATTERN_FH, ">$output_file") || die "Could not open $output_file for writing pattern file.\n";

	my $id=0;
	foreach my $oligo(@{$oligo_arr_ref}){
		print PATTERN_FH ">$PATPREFIX$id <mismatch=$error_rate>\n";
		print PATTERN_FH "$oligo\n";
		$id++;
	}
}

###############################################################################

sub load_oligos{
	# Loads sequences from file into hash
	my $filename=shift;
	my %hash;
	open(FH, "<$filename") || die "Could not open $filename\n";
	while(<FH>){
		chomp;
		my ($primers, $id)=split /\t/, $_;
		if(!($primers=~/^>/)){
			$hash{$primers}=$id;
		}
	}
	return(\%hash);
}

###############################################################################

sub parse_fuzznuc{

	my $search_output=shift;
	my %results_hash;
	my %sequence_len_hash;

	open(FZNC_FH, "<$search_output") || die "Could not open $search_output\n";

	my $current_seq_id;
	my $seq_length;
	while(<FZNC_FH>){
		chomp;
		if($_=~/^# Sequence: (\S+)\s+from:\s+(\d+)\s+to:\s+(\d+)/){
			$current_seq_id=$1;
			$seq_length=$3;

			# Placeholder in case there are no hits
			@{$results_hash{$current_seq_id}}=();	

			# I don't like reading this here, but it'll be faster than parsing the fasta file
			$sequence_len_hash{$current_seq_id}=$seq_length;	

		}elsif($_=~/^#/ || $_=~/^$/){
			# Noop
		}elsif($_=~/^  Start\s+End\s+Pattern_name\s+Mismatch\s+Sequence/){
			# Noop
		}else{
			$_=~s/^\s+//;
			my ($start, $end, $oligo_name, $mismatch_count, $sequence)=split /\s+/, $_;

			my $ori=1;
			if($start>$end){
				($start, $end)=($end, $start);
				$ori=-1;
			}

			$start--; #convert to 0-space based coordinates 	

			$mismatch_count=($mismatch_count eq ".")?0:int($mismatch_count);

			if(!($oligo_name=~s/^$PATPREFIX//)){
				die "Error Parsing oligo name out of fuzznuc for $_\n";
			}

			my $record=join ";", ($start, $end, $ori, $oligo_name, $mismatch_count);
			push @{$results_hash{$current_seq_id}}, $record;
		}

	}

	close(FZNC_FH);
	return(\%results_hash, \%sequence_len_hash);
}

###############################################################################

sub load_reads{
	my $input_fasta_file=shift;
	my %read_hash;
	
	open(FH, "<$input_fasta_file") || die "Could not read fasta file $input_fasta_file\n";

	while(<FH>){
		if($_=~/^>(\S+)/){
			$read_hash{$1}=1;
		}
	}
	return(\%read_hash);
}
