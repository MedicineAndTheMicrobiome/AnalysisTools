#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin ();
use lib "$FindBin::Bin";
use Getopt::Std;
use FileHandle;
use vars qw($opt_l $opt_p $opt_i $opt_r $opt_o);
use Sys::Hostname;
use File::Basename;

my $MOTHUR_BIN=`which mothur`;
chomp $MOTHUR_BIN;
my $MOTHUR_LOG="log";
my $TMP=hostname() . "." . $$;
my $NUM_PROC=8;


my $MATE_SCREEN_BIN="$FindBin::Bin/pipeline_utilities/Screen_MergedMates.r";
my $EXTRACT_SEQ_BIN="$FindBin::Bin/pipeline_utilities/Extract_Record_From_FASTA_By_List.pl";

getopts("l:pio:r:");
my $usage = "usage: 

$0 

	-l <fastq list filename>
	[-p (for a pair of fastq files)]
	[-i (for an interleaved fastq file)]
	[-r <relative path to append to fastq file lists>]
	-o <output directory>

	This script will call Mothur's make.contig() function
	on all the pairs of files specified.  If the input
	file is already paired, the fastq file will be split
	before joining the pairs.

	For the -p option, two columns should be specified.
	For the -i option, a single column should be specified.

	In both cases, the columns should contain the path of the fastq file.

	Using mothur from: $MOTHUR_BIN
	Using tmp root: $TMP
	Using mate screen from: $MATE_SCREEN_BIN

";

if(!(
	(defined($opt_p) || defined($opt_i))&& 
	defined($opt_o))){
	die $usage;
}

my $pairs_fn;
my $num_col;

if(defined($opt_p)){
	$pairs_fn=1;
	$num_col=2;
}else{
	$pairs_fn=0;
	$num_col=1;
}

my $path;
if(defined($opt_r)){
	$path=$opt_r;
}else{
	$path=".";
}

my $output_dir=$opt_o;
my $file_list=$opt_l;

if(!(-e $output_dir)){
	mkdir $output_dir;
}
if(!(-e $output_dir)){
	die "Could not make or find $output_dir.\n";
}

print STDERR "File List: $file_list\n";
print STDERR "Output Dir: $output_dir\n";
print STDERR "Paired? $pairs_fn\n";
print STDERR "Using temp file name root: $TMP\n";
print STDERR "Num Processors: $NUM_PROC\n";

###############################################################################

sub load_file_list{
	my $fn=shift;
	my $num_col=shift;

	my @files;
	open(FH, "<$fn") || die "Could not open $fn\n";
	while(<FH>){
		chomp;
		my @cols=split "\t", $_;
	
		if($num_col==1){
			push @files, $cols[0];
		}elsif($num_col==2){
			push @files, "$cols[0]\t$cols[1]";
		}
	}
	close(FH);
	return(\@files);
}

sub run_make_contigs{
	my $forward=shift;
	my $reverse=shift;
	my $output_dir=shift;

	my $exec_string=
		"$MOTHUR_BIN \"#set.logfile(name=$output_dir/$MOTHUR_LOG);make.contigs(ffastq=$forward, rfastq=$reverse, processors=$NUM_PROC);\"";
	#print STDERR "'$exec_string'\n";
	my $res=`$exec_string`;
	#print $res;

	if(!(-e $output_dir)){
		mkdir $output_dir;
	}
	if(!(-e "$output_dir/scrap")){
		mkdir "$output_dir/scrap";
	}
	if(!(-e "$output_dir/raw_merged_fasta")){
		mkdir "$output_dir/raw_merged_fasta";
	}
	if(!(-e "$output_dir/report")){
		mkdir "$output_dir/report";
	}
	if(!(-e "$output_dir/qual")){
		mkdir "$output_dir/qual";
	}
	if(!(-e "$output_dir/links")){
		mkdir "$output_dir/links";
	}
	if(!(-e "$output_dir/screen")){
		mkdir "$output_dir/screen";
	}
	if(!(-e "$output_dir/screened_fasta")){
		mkdir "$output_dir/screened_fasta";
	}


	my $new_root=$forward;
	$new_root=~s/\.fastq$//;
	my ($fname, $path)=File::Basename::fileparse($new_root);


	# Organize the output from the mate merge into separate directories
	`mv $new_root.scrap.contigs.* $output_dir/scrap`;
	`mv $new_root.trim.contigs.fasta $output_dir/raw_merged_fasta`;
	`mv $new_root.contigs.report $output_dir/report`;
	`mv $new_root.trim.contigs.qual $output_dir/qual`;
	`mv $forward $reverse $output_dir/links`;
	
	# Execute mate report analysis
	print STDERR "Performing merged mate screening...\n";
	$exec_string=
		"$MATE_SCREEN_BIN -i $output_dir/report/$fname.contigs.report -o $output_dir/screen/$fname.screen";
	print STDERR "'$exec_string'\n";
	$res=`$exec_string`;
	print STDERR "$res\n";

	# Extracted keep.list from raw merged
	print STDERR "Extracting merged mates that passed screening...\n";
	$exec_string=
		"$EXTRACT_SEQ_BIN " . 
		"-f $output_dir/raw_merged_fasta/$fname.trim.contigs.fasta " .
		"-l $output_dir/screen/$fname.screen.keep.list " .
		"| sed 's/\\t\$//' " .
		"> $output_dir/screened_fasta/$fname.fasta";
	print STDERR "'$exec_string'\n";
	$res=`$exec_string`;
	print STDERR "$res\n";

	

	# Makes:
	# s1_split.R1.contigs.report
	# s1_split.R1.scrap.contigs.fasta
	# s1_split.R1.scrap.contigs.qual
	# s1_split.R1.trim.contigs.fasta
	# s1_split.R1.trim.contigs.qual
	#
}

sub split_fastq{
	my $paired_fastq=shift;
	my $out_F=shift;
	my $out_R=shift;

	print STDERR "Splitting $paired_fastq...\n";
	if($paired_fastq=~/\.gz$/){
		open(FASTQ_FH, "zcat $paired_fastq| ") || die "Could not open $paired_fastq\n";
	}else{
		open(FASTQ_FH, "<$paired_fastq") || die "Could not open $paired_fastq\n";
	}

	open(OUT_F, ">$out_F") || die "Could not open $out_F\n";
	open(OUT_R, ">$out_R") || die "Could not open $out_R\n";

	my ($r1_id, $r1_seq, $r1_plus, $r1_qv);
	my ($r2_id, $r2_seq, $r2_plus, $r2_qv);

	my $num_recs=0;
	while(!eof(FASTQ_FH)){
		if($num_recs%2){
			$r2_id=<FASTQ_FH>;
			$r2_seq=<FASTQ_FH>;
			$r2_plus=<FASTQ_FH>;
			$r2_qv=<FASTQ_FH>;

			print OUT_F "$r1_id";
			print OUT_F "$r1_seq";
			print OUT_F "$r1_plus";
			print OUT_F "$r1_qv";

			print OUT_R "$r2_id";
			print OUT_R "$r2_seq";
			print OUT_R "$r2_plus";
			print OUT_R "$r2_qv";
		}else{
	               	$r1_id=<FASTQ_FH>;
			$r1_seq=<FASTQ_FH>;
			$r1_plus=<FASTQ_FH>;
			$r1_qv=<FASTQ_FH>;
		}	

		$num_recs++;
	}

	if($num_recs%2){
		die "Error:  Odd number of FASTQ records: $num_recs.\n";
	}

	close(FASTQ_FH);

	close(OUT_F);
	close(OUT_R);
}

###############################################################################
if($pairs_fn){
	print STDERR "Performing Analysis on F and R fastq files...\n";
}else{
	print STDERR "Performing Analysis on interleaved fastq file.\n";
}
print STDERR "\n";

# Read in file list
my $file_list_ref=load_file_list($file_list, $num_col);

# For each file run make_contigs in Mothur
foreach my $files(@{$file_list_ref}){
	my ($for_fname, $rev_fname);

	print STDERR "$files\n";
	
	my ($fname, $rname);

	if(!$pairs_fn){
		my ($unsplit_fname, $fpath)=fileparse($files);
		$unsplit_fname=~s/\.fastq$//;

		# Append path if relative
		if(substr($files,0,1) ne "/"){
			$files="$path/$files";
		}

		# Create temp file names for split results
		my $temp_for_fname="$output_dir/tmp.$TMP.for.fastq";	
		my $temp_rev_fname="$output_dir/tmp.$TMP.rev.fastq";	

		# Split interleaved file
		split_fastq($files, $temp_for_fname, $temp_rev_fname);

		# Assign name to split file
		$fname="$unsplit_fname.for.fastq";
		$rname="$unsplit_fname.rev.fastq";

		# Link to assigned names to temp split files
		symlink $temp_for_fname, "$output_dir/$fname";
		symlink $temp_rev_fname, "$output_dir/$rname";
		
	}else{

		# Set up files if they are already split into F and R
		my ($for_file, $rev_file)=split "\t", $files;	

		# Append path if relative
		if(substr($for_file,0,1) ne "/"){
			$for_fname="$path/$for_file";
		}else{
			$for_fname=$for_file;
		}

		if(substr($rev_file,0,1) ne "/"){
			$rev_fname="$path/$rev_file";
		}else{
			$rev_fname="$rev_file";
		}

		# Get root file name
		my ($fpath, $rpath);
		($fname, $fpath)=fileparse($for_fname);
		($rname, $rpath)=fileparse($rev_fname);

		# Link to original files
		symlink $for_fname, "$output_dir/$fname";
		symlink $rev_fname, "$output_dir/$rname";

	}
	
	run_make_contigs("$output_dir/$fname", "$output_dir/$rname", $output_dir);
	
	print STDERR "ok...\n";
}

###############################################################################

print STDERR "done.\n";

