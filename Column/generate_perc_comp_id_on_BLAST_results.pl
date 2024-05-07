#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Basename;

###############################################################################

use vars qw ($opt_s $opt_q $opt_p $opt_t $opt_r $opt_o);
getopts("s:q:p:t:r:o:");

my $usage = "
	$0
	    [-s <Subject (reference) fasta file, used in blast>]
	    [-q <Query (reads) fasta file, used in blast>]

	    [-t <Subject (reference) length file>]
	    [-r <Query (reads) length file>]
	
	    [-o <output directory>]

	    -p \"<blast output file pattern>\"

	The output will go to <input blast output>.comp_id

	This program will read in the blast output when the -m 8 option 
	has been used, and then compute the composite percent identities
	for every aligment.  
	
	The composite id of query is: (alignment length)/(query length)*(alignment percent id)
	The composite id of sbjct is: (alignment length)/(sbjct length)*(alignment percent id)
		
	Order of the output:
		Query_id (read id)
		Comp_perc_of_query
		Query_length

		Subject_id (reference id)
		Comp_perc_of_subject
		Subject_length

		perc_identity
		alignment_length
		bit_score
		evalue


	The fasta files are just used for their length.  Make sure
	-s and -q don't have any different sequences with the same
	id, because they are read into the same sequence id to length
	hash.

	If you have the lengths already, just use the -t and -r options.

	if the -o <output directory> is specified then the output .comp_id
	file will be placed into the <output directory>.  Else the output
	will be placed in the same directory as the blast output (input).

";

if(!defined($opt_p)){
	die $usage;
}

my $subject_fasta=$opt_s;
my $query_fasta=$opt_q;
my $blastoutput_pat=$opt_p;
my $subject_lengths=$opt_t;
my $query_lengths=$opt_r;
my $comp_id_output_dir=$opt_o;

###############################################################################

my %seq_length_hash;
my $count=0;

###############################################################################

my ($defline, $prev_defline, $sequence);

if(defined($subject_fasta) || defined($query_fasta)){

	open(FASTA_FH, "cat $subject_fasta $query_fasta | ") || 
		die "Could not open one of the sequence fasta files\n";
	print STDERR "Processing FASTA file...\n";
	while(<FASTA_FH>){
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
	print STDERR "\n";
	close(FASTA_FH);

}

if(defined($subject_lengths) || defined($query_lengths)){

	open(LENGTHS_FH, "cat $subject_lengths $query_lengths | ") || 
		die "Could not open one of the length files\n";

	while(<LENGTHS_FH>){
		chomp;
		my ($id, $length)=split /\t/, $_;
		if($id=~/^>/){
			$id=~s/^>//;
		}
		$seq_length_hash{$id}=$length;	
	}

	close(LENGTHS_FH);
}


sub process_record{
	my $defline=shift;
	my $sequence=shift;

	$count++;
	if(!($count %10000)){
		print STDERR ".";
	}

	my $id;
	my $len=length($sequence);

	if($defline=~/^>(\S+)/){
		$id=$1;
		$seq_length_hash{$id}=$len;
	}else{
		die "Could not parse out ID, $_\n"
	}
}


###############################################################################

my @files=split /\n/, `ls -1 $blastoutput_pat`;

foreach my $file(@files){

	print STDERR "Processing $file\n";

	# If file is compressed, use zcat, else just redirect in
	if($file=~/\.gz$/){
		open(FH, "zcat $file |") || die "Could not open '$file' (compressed)\n";
		$file=~s/\.gz$//;
	}else{
		open(FH, "<$file") || die "Could not open '$file' (uncompressed)\n";
	}

	# If output is defined, move output to specified output dir
	if(defined($comp_id_output_dir)){
		my ($name, $path)=fileparse($file);
		$file="$comp_id_output_dir/$name";
	}
	
	print STDERR "Writing to: $file\.comp_id\n";
	open(OUT, ">$file\.comp_id") || die "Could not open '$file\.comp_id'\n";

	while(<FH>){
		chomp;

		if($_=~/^#/){
			# Skip comments.
			next;
		}

	 	my ($Query_id,$Subject_id,$perc_identity,
			$alignment_length,$mismatches,$gap_openings,
			$qstart,$qend,$sstart,$send,$evalue,$bit_score)=split /\t/, $_;

		my $subject_comp_id="NA";
		my $query_comp_id="NA";

		my $sub_seq_length="NA";
		my $qry_seq_length="NA";

		if(!defined($seq_length_hash{$Subject_id})){
			$sub_seq_length="NA";
			$subject_comp_id="NA";
		}else{
			$sub_seq_length=$seq_length_hash{$Subject_id};
			#$subject_comp_id=sprintf("%3.2f", $perc_identity*$alignment_length/$sub_seq_length);
			$subject_comp_id=sprintf("%3.2f", 
				$perc_identity*(abs($send-$sstart)+1)/$sub_seq_length);
		}

		if(!defined($seq_length_hash{$Query_id})){
			$qry_seq_length="NA";
			$query_comp_id="NA";
		}else{
			$qry_seq_length=$seq_length_hash{$Query_id};
			#$query_comp_id=sprintf("%3.2f", $perc_identity*$alignment_length/$qry_seq_length);
			$query_comp_id=sprintf("%3.2f", 
				$perc_identity*(abs($qend-$qstart)+1)/$qry_seq_length);
		}
		

		my $outstr=join "\t", (
			$Query_id,
			$query_comp_id,
			$qry_seq_length,

			$Subject_id,
			$subject_comp_id,
			$sub_seq_length,

			$perc_identity,
			$alignment_length,
			$bit_score,
			$evalue			
		);

		print OUT "$outstr\n";
	}
}

