#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin ();
use Getopt::Std;
use FileHandle;
use File::Basename;
use Sys::Hostname;

my $PIPELINE_UTIL_PATH="$FindBin::Bin/../../pipeline_utilities";
my $EXTRACT_FASTA_PATH="$PIPELINE_UTIL_PATH/Extract_Record_From_FASTA_By_List.pl";

print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";

my $OTU_CUTOFF=0.03;
my $NUM_HITS=2000;
my $TMP_PATH="/usr/local/scratch";

use vars qw($opt_f $opt_r $opt_o $opt_m $opt_c $opt_s $opt_t);
getopts("f:r:o:m:c:s:t:");
my $usage = "usage: 

$0 

	-f <fasta file>
	-r <reference fasta, formatdb'd and curated>
	-o <output root>
	-m <accession to species map>
	[-c <similarity cutoff, def=$OTU_CUTOFF>]
	[-s <sample id>]
	[-t <tmp path, def=$TMP_PATH>]

	This script will take in a query fasta file,
	presumably representing the reads from a single
	OTU and try to estimate the proportions of
	the species underlying it based on the
	reference fastas file.

	For each read, if no references are within the cutoff 
	from it then no subject will be associated with 
	it.  If there is more than 1 best hit and there is
	no consensus species name associated with it
	then it will be marked ambiguous.

	The Output file will be a single line summary_table:

	sample_id total species1 species2 ... speciesN

";

if(!(
	defined($opt_f) && 
	defined($opt_r) && 
	defined($opt_o))){
	die $usage;
}

my $qry_fasta=$opt_f;
my $blast_db=$opt_r;
my $cutoff=$OTU_CUTOFF;
my $output_root=$opt_o;

if(defined($opt_c)){
	$cutoff=$opt_c;
}

my $sample_id=$opt_s;
if(!defined($opt_s)){
	my ($name, $path, $ext)=fileparse($qry_fasta, qr/\.[^.]*/);
	$sample_id=$name;
}

my $map_file="";
if(defined($opt_m)){
	$map_file=$opt_m;
}

my $tmppath=$TMP_PATH;
if(defined($opt_t)){
	$tmppath=$opt_t;
}
my $host_id=hostname;
my $tmp_fn="$tmppath/$host_id\.$$.blastout";



print STDERR "Query FASTA: $qry_fasta\n";
print STDERR "Blast DB Reference: $blast_db\n";
print STDERR "Map File: $map_file\n";
print STDERR "Cutoff: $cutoff\n";
print STDERR "Sample ID: $sample_id\n";
print STDERR "Output Root: $output_root\n";
print STDERR "Temp File: $tmp_fn\n";
print STDERR "\n";
print STDERR "Num Hits to Analyze: $NUM_HITS\n";

###############################################################################

###############################################################################

sub get_seq_lengths{
	my $fasta=shift;
	my %length_hash;

	sub process_record{
		my $defline = shift;
		my $sequence = shift;

		my $length=length($sequence);
		if($defline=~/^>(\S+)/){
			$length_hash{$1}=$length;
		}else{
			die "Error: processing defline for sequence id.\n";	
		}
	}


	open(FH, "<$fasta") || die "Could not open $fasta for sequence length counting.\n";
	my ($defline, $prev_defline, $sequence);
	while(<FH>){
		chomp;

		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				process_record($prev_defline, $sequence);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	process_record($prev_defline, $sequence);
	close(FH);

	my $num_lengths_loaded=scalar keys %length_hash;
	print STDERR "Number of Lengths Loaded: $num_lengths_loaded\n";

	return(\%length_hash);
}

###############################################################################

sub load_seq_map{
	my $fn=shift;
	my %hash;

	open(FH, "<$fn") || die "Could not open $fn\n";
	while(<FH>){
		chomp;
		my ($src, $dst)=split "\t", $_;
		$hash{$src}=$dst;
	}
	close(FH);

	return(\%hash);
}

###############################################################################

sub blast_sequences{
	my $qry_fasta=shift;
	my $database_fasta=shift;
	my $tmp_blast_out=shift;
	my $num_hits=shift;

	my $cmd="blastall " .
		"-p blastn " . 
		"-d $database_fasta " . 
		"-i $qry_fasta " .
		"-e 1e-100 " .
		"-m 8 " .
		"-b $num_hits " .
		"-o $tmp_blast_out " .
		"-F F ";

	print STDERR "Executing: $cmd\n";
	system($cmd);
	print STDERR "ok.\n";

}

###############################################################################

#1 #Query
#2 #Subject
#3 #% id
#4 #alignment length
#5 #mistmatches
#6 #gap openings
#7 #q.start
#8 #q.end
#9 #s.start
#10 #s.end
#11 #e-value
#12 #bit score

my $QRY_ID=0;
my $SUBJ_ID=1;
my $PERC_ID=2;
my $ALN_LEN=3;
my $BIT_SCR=11;

sub classify_sequences{
	my $blast_out_fn=shift;
	my $len_hash_ref=shift;
	my $taxa_hash_ref=shift;
	my $cutoff=shift;
	my $max_recs=shift;

	my $sim_thresh=(1.0-$cutoff)*100.0;

	sub process_blast_record{
		my $arr_ref=shift;
		my $qry_id=${${$arr_ref}[0]}[0];

		# All records here should have the same query ID

		# Find highest bitscore
		my $max_bitscr=0;
		foreach my $rec_ref(@{$arr_ref}){
			my $cur_bitscr=${$rec_ref}[$BIT_SCR];
			if($max_bitscr<$cur_bitscr){
				$max_bitscr=$cur_bitscr;
			} 
		}
		
		# Only keep highest bitscore alignments with close percent similarity
		my @best_records=();
		my %subjects_hit_hash;
		foreach my $rec_ref(@{$arr_ref}){
                	my $cur_bitscr=${$rec_ref}[$BIT_SCR];
			if($cur_bitscr == $max_bitscr){

				push @best_records, $rec_ref;

				my $subj_id=${$rec_ref}[$SUBJ_ID];
				my $qry_id=${$rec_ref}[$QRY_ID];
				my $perc_id=${$rec_ref}[$PERC_ID];
				my $aln_len=${$rec_ref}[$ALN_LEN];

				# Compute percent composite identity
				my $qry_len=${$len_hash_ref}{$qry_id};
				if($aln_len>$qry_len){
					$aln_len=$qry_len;
				}
				my $perc_comp_id=$perc_id*($aln_len/$qry_len);

				# Only keep records close to the the reference
				my $taxa=${$taxa_hash_ref}{$subj_id};
				if($perc_comp_id >= $sim_thresh){
					if(!defined($subjects_hit_hash{$taxa})){
						$subjects_hit_hash{$taxa}=1;
						print STDERR "$subj_id [$perc_comp_id%]\n";
					}else{
						$subjects_hit_hash{$taxa}++;
					}
				}
			}
		}

		# Make sure the number of best ties is less than the number of records.
		my $num_recs=$#best_records+1;
		if($num_recs==$max_recs){
			print STDERR "WARNING: Max records hit: $max_recs.\n";
			print STDERR "Either reduce database size or increase num blast records returned.\n";
		}		

		my @hit_arr;
		foreach my $subj(keys %subjects_hit_hash){
			push @hit_arr, "$subj";				
			#push @hit_arr, "$subj($subjects_hit_hash{$subj})";				
		}

		# Determine if there were any hits
		my $hit_str;
		if($#hit_arr>-1){
			$hit_str=join "/", @hit_arr;		
		}else{
			$hit_str="NoCloseRef";
		}
		print STDERR "$qry_id:\t$hit_str\n";

		# Return single string analysis
		return("$qry_id\t$hit_str");

	}


	my @class_results;
	
	# Load blast output, one group at a time.
	open(FH, "<$blast_out_fn") || die "Could not open $blast_out_fn for analysis.\n";
	my ($prev_id, $cur_id)=("","");
	my @cur_records=();
	my $taxa;	
	while(<FH>){
		chomp;

		my @cols=split "\t", $_;
		
		my $cur_id=$cols[$QRY_ID];
		if(($cur_id ne $prev_id) && ($prev_id ne "")){
			$taxa=process_blast_record(\@cur_records);
			push @class_results, $taxa;
			@cur_records=();
		}	
		$prev_id=$cur_id;

		push @cur_records, \@cols;
	}
	$taxa=process_blast_record(\@cur_records);
	push @class_results, $taxa;
	close(FH);

	return(\@class_results);
}

###############################################################################

sub output_summary_table{
	my $output_root=shift;
	my $class_arr_ref=shift;
	my $sample_id=shift;
	my $total_seq=shift;

	my %category_hash;

	open(FH_ST, ">$output_root.summary_table.tsv") || die "Could not open $output_root.summary_table.tsv\n";

	my $subtotal=0;
	foreach my $class(@{$class_arr_ref}){
		my ($qry_id, $taxa)=split "\t", $class;	
		if(!defined($category_hash{$taxa})){
			$category_hash{$taxa}=1;
		}else{
			$category_hash{$taxa}++;
		}
		$subtotal++;
	}
	
	my @categories=sort keys %category_hash;
	my $hdr=join "\t", ("sample_id", "total", @categories, "NoReference");
	print FH_ST "$hdr\n";
	print FH_ST "$sample_id\t$total_seq";

	foreach my $cat(@categories){
		print FH_ST "\t$category_hash{$cat}";
	}
	print FH_ST "\t" . ($total_seq-$subtotal);

	print FH_ST "\n";

	close(FH_ST);
	return;
}


sub output_classification_detail{
	my $output_root=shift;
	my $class_arr_ref=shift;

	open(FH_CD, ">$output_root.classify.detail.tsv") || die "Could not open $output_root.classify.detail.tsv\n";

	foreach my $class(@{$class_arr_ref}){
		print FH_CD "$class\n";
	}

	close(FH_CD);
	return;
}

###############################################################################

print STDERR "Loading sequence lengths...\n";
my $seq_len_hash_ref=get_seq_lengths($qry_fasta);

print STDERR "Loading taxa map...\n";
my $seq_to_taxa_hash_ref=load_seq_map($map_file);

print STDERR "Performing blast...\n";
blast_sequences($qry_fasta, $blast_db, $tmp_fn, $NUM_HITS);

#$tmp_fn="/usr/local/scratch/METAGENOMICS/kli/tmp/lserver2.30783.blastout";

print STDERR "Classifying sequences...\n";
my $class_arr_ref=classify_sequences(
	$tmp_fn,
	$seq_len_hash_ref,
	$seq_to_taxa_hash_ref,
	$cutoff,
	$NUM_HITS
);

my $num_total_seq=scalar keys %{$seq_len_hash_ref};
output_summary_table($output_root, $class_arr_ref, $sample_id, $num_total_seq);
output_classification_detail($output_root, $class_arr_ref);

system("rm $tmp_fn");

###############################################################################

print STDERR "done.\n";

