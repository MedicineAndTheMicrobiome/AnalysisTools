#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Std;
use vars qw($opt_d $opt_s $opt_p $opt_o $opt_i $opt_t);
use File::Basename;

my $RAW_FASTQ_FILES="Run";

getopts("d:s:p:o:i:t");
my $usage = "usage: 

$0 

	Required:
	-d <SequencingRuns Directory>
	-s <'list of sequencing runs to search' file>
	-p <'list of project IDs to include', 0000 and 0001, negative controls automatically included>
	-o <Output directory>

	Optional:
	[-i <include list, exclude other samples found>]
	[-t <build tgz file of fastq files, otherwise just create directory of symbolic links>]

	This script will look into the SequencingRuns directory, finding the sequencing runs
	specified in the sequencing runs list.  For each sequencing run, it will look into
	
	$RAW_FASTQ_FILES

	The paired fastq files will be copied into the output directory.

	A single tsv file will be generated which tries to guess the above variables
	based on the sample ID.
	
";

if(
	!defined($opt_d)||
	!defined($opt_s)||
	!defined($opt_p)||
	!defined($opt_o)
){
	die $usage;
}

my $SequencingRunsDirectory=$opt_d;
my $SequencingRunsListFname=$opt_s;
my $ProjectIDsListFname=$opt_p;
my $OutputDirectory=$opt_o;

my $IncludeList="";
if(defined($opt_i)){
	$IncludeList=$opt_i;
}

my $BuildTGZ=0;
if(defined($opt_t)){
	$BuildTGZ=1;
}

###############################################################################

print STDERR "\n";
print STDERR "Sequencing Runs Directory: $SequencingRunsDirectory\n";
print STDERR "Sequencing Run Target List File Name: $SequencingRunsListFname\n";
print STDERR "Project IDs List File Name: $ProjectIDsListFname\n";
print STDERR "Output Directory: $OutputDirectory\n";
print STDERR "\n";
print STDERR "Include List: $IncludeList\n";
print STDERR "Build TGZ: $BuildTGZ\n";

###############################################################################

sub load_file_list{
	my $list=shift;

	system("dos2unix $list");

	open(IN_FH, "<$list") || die "Could not open $list\n";
	my @load_file;
	while(<IN_FH>){
		chomp;
		if(substr($_, 0, 1) eq "#"){ next;} # Skip comments
		push @load_file, $_;
	}
	close(IN_FH);	
	return(\@load_file);
}

###############################################################################

sub make_dir{
	my $dir=shift;
	print STDERR "Making $dir.\n";
	if(!(-e $dir)){
		mkdir $dir;
	}else{
		print STDERR "$dir exists.\n";
	}
	if(!(-e $dir)){
		die "Could not make or find $dir\n";
	}
}

###############################################################################
# Load targeted sequence runs

my $target_seq_runs_ref=load_file_list($SequencingRunsListFname);
my $target_projects_ids_ref=load_file_list($ProjectIDsListFname);

unshift @{$target_projects_ids_ref}, ("0000", "0001", "0002");


# Clean up project ids a little and make sure they are unique
my %pid_hash;
foreach my $pid(@{$target_projects_ids_ref}){
	if(length($pid)>1){
		$pid_hash{$pid}=1;
	}
}
@{$target_projects_ids_ref}=sort keys %pid_hash;

###############################################################################

print STDOUT "\nIdentify Existence of Targeted Sequencing Runs:\n\n";

my $missing_sr_dirs=0;
foreach my $sr(@{$target_seq_runs_ref}){
	#print STDOUT "\t$sr\n";

	if(-e "$SequencingRunsDirectory/$sr"){
		print "$sr found.\n";
	}else{
		print "$sr NOT FOUND!!!\n";
		$missing_sr_dirs=1;
	}
}

if($missing_sr_dirs){
	print STDERR "\n";
	print STDERR "ERROR: Missing sequencing run directories.\n";
	print STDERR "Please correct error(s).\n";
	print STDERR "(You man need to fix typos in the sequencing\n";
	print STDERR " names, or remove them from the list to proceed.)\n";
	print STDERR "\n";
}

###############################################################################

my $dont_copy=0;


# Try to make output dir
make_dir($OutputDirectory);
make_dir("$OutputDirectory/fastq");

my @run_list;
my @file_list;
my @f1_list;
my @f2_list;

foreach my $sr(@{$target_seq_runs_ref}){
	
	my $seqrun_full_dir="$SequencingRunsDirectory/$sr";

	print STDERR "Looking for target files in: $SequencingRunsDirectory/$sr\n";

	foreach my $projid (@{$target_projects_ids_ref}){

		print STDERR "Looking for project ID: $projid\n";

		my @files=`ls $seqrun_full_dir/$RAW_FASTQ_FILES/$projid-*_R1_*.fastq.gz`;

		foreach my $file(@files){
			chomp $file;
			print STDERR "Copying file:\n\t$file\nto\n\t$OutputDirectory/$sr\n\n";

			my $f1=$file;
			my $f2=$f1;
			$f2=~s/_R1_001.fastq.gz/_R2_001.fastq.gz/;

			my $f1_fname=basename($f1);
			my $f2_fname=basename($f2);

			print STDERR "$f1_fname\n";
			print STDERR "$f2_fname\n";

			my $f1_newname="$sr\_$f1_fname";
			my $f2_newname="$sr\_$f2_fname";

			print STDERR "$f1_newname\n";
			print STDERR "$f2_newname\n";

			push @run_list, $sr;
			push @file_list, $f1_fname;
			push @f1_list, $f1_newname;
			push @f2_list, $f2_newname;

			if(!$dont_copy){
				print STDERR `cp -s $f1 $OutputDirectory/fastq/$f1_newname`;
				print STDERR `cp -s $f2 $OutputDirectory/fastq/$f2_newname`;
			}
		}

	}

}

my $num_entries=$#run_list+1;
print STDERR "Num Entries: $num_entries\n";

my %sample_type_hash;

$sample_type_hash{"ST"}="Stool#Human";
$sample_type_hash{"SINUS"}="Sinus#Human";
$sample_type_hash{"OW"}="OralWash#Human";
$sample_type_hash{"SP"}="Sputum#Human";
$sample_type_hash{"NP"}="NasopharyngealSwab#Human";
$sample_type_hash{"BAL"}="BronchoalveolarLavage#Human";
$sample_type_hash{"SKIN"}="Skin#Human";
$sample_type_hash{"SAL"}="Saliva#Human";
$sample_type_hash{"MST"}="Stool#Mouse";
$sample_type_hash{"MLUNG"}="Lung#Mouse";
$sample_type_hash{"RST"}="Stool#Rat";

my %medium_hash;
$medium_hash{"Stool"}="feces [UBERON:0001988]";
$medium_hash{"OralWash"}="saliva [UBERON:0000165]";
$medium_hash{"Saliva"}="saliva [UBERON:0000165]";

my %local_scale_hash;
$local_scale_hash{"Stool"}="colon [UBERON:0001155]";
$local_scale_hash{"OralWash"}="mouth [UBERON:0001165]";
$local_scale_hash{"Saliva"}="mouth [UBERON:0001165]";

my %control_hash;
$control_hash{"0000"}="Extraction_Kit_Negative_Control#kit_reagents#metagenome";
$control_hash{"0001"}="PCR_Amplification_Negative_Control#pcr_reagents#metagenome";
$control_hash{"0002"}="PCR_Amplification_PositiveZymo_Control#control_DNA#synthetic metagenome";

my %organism_hash;
$organism_hash{"Human"}="Homo sapiens";
$organism_hash{"Rat"}="Rattus norvegicus";
$organism_hash{"Mouse"}="Mus musculus";

my $env_broad_scale_default="multicellular organism [UBERON:0000468]";

my @header=(
	"sample_name",
	"library_id",
	"title",
	"organism",
	"collection_date", 
	"env_broad_scale",
	"env_local_scale",
	"env_medium",
	"host",
	"host_subject_id1",
	"host_subject_id2",
	"filename",
	"filename2"
);


my $outfile="$OutputDirectory/SRA_fields.tsv";
open(FH, ">$outfile") || die "Could not open $outfile.\n";

print FH (join "\t", @header) . "\n";

for(my $i=0; $i<$num_entries; $i++){
	print STDERR "\n\n";
	print STDERR "INPUT: $run_list[$i],$file_list[$i]\t$f1_list[$i]/$f2_list[$i]\n";

	# Parse file name
	my @file_toks=split "_", $file_list[$i];
	my $sample_name=$file_toks[0];
	my @sample_name_toks=split "-", $sample_name;

	# Parse Dates
	my $sample_date_raw=$sample_name_toks[3];
	my $sample_year=substr($sample_date_raw, 0, 4);
	my $sample_mon=substr($sample_date_raw, 4, 2);
	my $sample_day=substr($sample_date_raw, 6, 2);
	my $sample_date_form="$sample_year-$sample_mon-$sample_day";

	# Parse sample types
	my $sample_type_rec=$sample_type_hash{$sample_name_toks[4]};
	my @sample_type_tok=split "#", $sample_type_rec;
	my $tissue=$sample_type_tok[0];
	my $env_medium=$medium_hash{$tissue};
	my $env_local_scale=$local_scale_hash{$tissue};
	my $env_broad_scale=$env_broad_scale_default;
	my $title="16S_Microbiome_" . $tissue;
	my $host=$organism_hash{$sample_type_tok[1]};
	my $organism=lc("$sample_type_tok[1] $tissue microbiome");

	my $host_subject_id1=$sample_name_toks[1];
	my $host_subject_id2=$sample_name_toks[2];

	my $library_id="$run_list[$i]\.$sample_name";

	# Override variable values if it is a control
	my $control_info=$control_hash{$sample_name_toks[0]};
	if(defined($control_info)){
		my ($control_name, $control_scale, $control_organism)=split "#", $control_info;
		$sample_name=$control_name;
		$env_local_scale=$control_scale;
		$env_broad_scale=$control_scale;
		$env_medium=$control_scale;
		$organism=$control_organism;
		$host=$control_scale;
		
		# Reformat control processing date
		$sample_date_raw=$sample_name_toks[2];
		$sample_year=substr($sample_date_raw, 0, 4);
		$sample_mon=substr($sample_date_raw, 4, 2);
		$sample_day=substr($sample_date_raw, 6, 2);
		$sample_date_form="$sample_year-$sample_mon-$sample_day";
	}

	# Generate output string
	my $outstr=join "\t", (
		$sample_name, 
		$library_id, 
		$title, 
		$organism,
		$sample_date_form, 
		$env_broad_scale,
		$env_local_scale, 
		$env_medium, 
		$host,	#scientific name
		$host_subject_id1, 
		$host_subject_id2, 
		$f1_list[$i], 
		$f2_list[$i]
	);

	print FH "$outstr\n";
	
}

close(FH);

###############################################################################

if($BuildTGZ){
	print STDERR "\n\n";
	print STDERR "Building tgz file...  This may take a while!\n";
	# the h option deferences links
	system("tar -cvzhf $OutputDirectory/fastq.tgz $OutputDirectory/fastq");
}else{
	print STDERR "\n\n";
	print STDERR "********************************************************\n";
	print STDERR "The option to build a fastq.tgz was not exercised.\n";
	print STDERR "Make sure you specify the -t option make it so!.\n";
	print STDERR "********************************************************\n";
	print STDERR "\n\n";
}	


###############################################################################

print STDERR "Done.\n";

