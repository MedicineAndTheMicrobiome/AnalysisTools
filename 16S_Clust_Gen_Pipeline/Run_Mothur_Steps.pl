#!/usr/bin/env perl

###############################################################################

use FindBin ();
use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_f $opt_g $opt_r $opt_o $opt_p $opt_c);

my $MOTHUR_BIN="/usr/bin/mothur";

my $PIPELINE_UTIL_PATH="$FindBin::Bin/pipeline_utilities";
print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";

my $OTU_TO_ST_BIN="$PIPELINE_UTIL_PATH/OTU_To_SummaryTable/Convert_MothurSharedFile_to_SummaryTable.r";
my $TAXA_TO_ST_BIN="$PIPELINE_UTIL_PATH/Taxonomy_To_SummaryTable/Convert_MothurTaxonomy_to_SummaryTable.pl";
my $COUNT_NAMES_BIN="$PIPELINE_UTIL_PATH/Count_Names/Count_Names.pl";

#my $CURRENT_16S_ALIGNMENT=
#	"/usr/local/devel/DAS/users/kli/SVN/DAS/16sDataAnalysis/trunk/16S_OTU_Generation/silva.nr_v119.align";

# Execution settings
my $TIMING_LOGNAME="timing_log.tsv";
my $MOTHUR_LOG="mothur.current.logfile";
my $COUNTS_LOGNAME="counts.logfile";
my $DEF_NPROC=4;
my $DEF_CLUST_CUTOFF=0.10;

###############################################################################

getopts("f:g:r:o:p:c:");
my $usage = "usage: 

$0 
	-f <16S fasta file, quality trimmed>
	-g <groups file, read-to-sample id file>
	-r <reference 16S alignments, e.g. [abs path]/silva.nr_v119.align >
	-o <output directory>

	[-p <num processors, default=$DEF_NPROC>]
	[-c <maximum cluster cutoff, default=$DEF_CLUST_CUTOFF>]

	This script will run through all the necessary steps to go from
	cleaned up FASTA sequences per sample, to OTU generation and taxonomic
	assignment.
	
	The (-f) 16S fasta file should quality, primer, adapter trimmed and length
	filtered.  This pipeline will not do any of that for you.  The fasta
	file should contain all the 16S sequence you want to cluster across
	all samples.

	The (-g) groups file contains a mapping from reads to sample id.
	You can generate this file with the script Assign_Sample_IDSs_To_Reads.pl.

	Make sure you specify the latest version of the 16S alignments for (-r).

	The (-o) output directory should be performed in scratch.  Then you
	can zip and archive the intermediate files.  You just want to make sure
	don't run this pipeline in a project directory where you may be charged
	for one time temporary disk usage.

	The number of processors can be as many as you have access to. It 
	will significantly improve the alignment to the reference step.


	You perform a test run based in the sample files in 'testing_files'.

	If perform a 'ls -ltr' in the output directory you will see the steps
	that have completed.  You will see a ##_<command> file for each 
	step in the pipeline.  These files are log files.  If you delete
	one of these log files, you can rerun the pipeline starting from that
	point automatically.

	You will also see:
		counts.logfile:	This contains a count of the reads and representatives
			that have been collapsed or removed during the pipeline
		timing_log.tsv: This contains a breakdown of the steps and how
			many CPU seconds and when the step was started.
		Summary_Tables: This contains the OTU and taxonomic summary tables
			across all samples.

";

if(!(
	defined($opt_f) && 
	defined($opt_g) && 
	defined($opt_r) && 
	defined($opt_o))){
	die $usage;
}

my @overall_begin_time=times;

my $input_fasta=$opt_f;
my $groups_file=$opt_g;
my $ref_16s_align=$opt_r;
my $output_dir=$opt_o;
my $num_proc=defined($opt_p)?$opt_p:$DEF_NPROC;
my $clust_cutoff=defined($opt_c)?$opt_c:$DEF_CLUST_CUTOFF;

print STDERR "Using Mothur at: $MOTHUR_BIN\n";
print STDERR "Input FASTA File: $input_fasta\n";
print STDERR "Groups File: $groups_file\n";
print STDERR "Reference 16S Alignments: $ref_16s_align\n";
print STDERR "Output Directory: $output_dir\n";
print STDERR "Num Processors: $num_proc\n";
print STDERR "Cluster Cutoff: $clust_cutoff\n";

if(!(-e $output_dir)){
	print STDERR "Making $output_dir...\n";
	mkdir $output_dir;
}else{
	print STDERR "$output_dir already exists...\n";
}
print "\n";

###############################################################################

sub get_elapsed_time{
	my $begin_ref=shift;
	my $end_ref=shift;

	my $begin_tot=0;
	my $end_tot=0;

	for(my $i=0; $i<4; $i++){
		$begin_tot+=${$begin_ref}[$i];
		$end_tot+=${$end_ref}[$i];
	}
	
	my $elapsed=$end_tot-$begin_tot;
	return($elapsed);
}

sub format_datetime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

	my $year_str=$year+1900;
	my $mon_str=$mon+1;
	my $day_str=$mday;
	my $date="$mon_str/$day_str/$year_str";
	
	my $time=sprintf("%02i:%02i:%02i", $hour,$min,$sec);
	return($date, $time);
}


my $timing_log="$output_dir/$TIMING_LOGNAME";

sub log_time{
	my $date=shift;
	my $time=shift;
	my $begin_time_ref=shift;
	my $end_time_ref=shift;
	my $step=shift;
	my $notes=shift;

	my $run_time=get_elapsed_time($begin_time_ref, $end_time_ref);	

	my $print_hdr=0;
	if(!-e $timing_log){
		$print_hdr=1;
	}

	open(LOG, ">>$timing_log") || die "Could not open $timing_log for appending.\n";
	if($print_hdr){
		print LOG "Start_Date\tStart_Time\tCommand\tCPU_Time\tNotes\n";
	}
	my $run_time_str=sprintf("%3.2f", $run_time);
	print LOG "$date\t$time\t$step\t$run_time_str\t$notes";
	print LOG "\n";
	close(LOG);

	return;
}

###############################################################################

sub check_logfile_for_errors{
	my $log_fname=shift;
	my $errors_found=0;
	open(FH, "<$log_fname") || die "Could not open $log_fname to check for errors.\n";

	while(<FH>){
		my $line=$_;
		if($line=~/\[ERROR\]/){
			$errors_found=1;
			last;
		}
	}

	close(FH);
	return($errors_found);
}


my $no_more_aborts=0;
my $step=1;

sub execute_mothur_cmd{
	my $cmd=shift;
	my $param=shift;

	$param=~s/\s+//g;	

	print STDERR "\n";
	print STDERR "***************************************************************\n";
	print STDERR "*  Executing $cmd in Mothur...\n";
	print STDERR "***************************************************************\n";

	my $step_str=sprintf("%02i", $step);
	my $logfile="$output_dir/$step_str\_$cmd";
	my $run_time=0;
	my $step_skipped=0;

	# Get start time
	my ($date, $time)=format_datetime();
	my $notes="";

	my (@begin_time, @end_time);
	if(!(-e $logfile) || $no_more_aborts){
	
		my $exec_string="$MOTHUR_BIN \"#set.logfile(name=$output_dir/$MOTHUR_LOG,append=T);$cmd($param)\"";
		print STDERR "$exec_string\n";

		@begin_time=times;
		my $out=`$exec_string > $logfile 2>&1`;
		@end_time=times;

		#open(FH, ">$logfile") || die "Could not open $logfile.\n";
		#print FH "$out\n";
		#close(FH);

		# Check for errors in the log file
		my $has_errors=check_logfile_for_errors($logfile);
		if($has_errors){
			print STDERR "\n\n";
			print STDERR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
			print STDERR "!                                                              !\n";
			print STDERR "!  TERMINATING:  Error found in mothur log file!!!             !\n";
			print STDERR "!                Could not complete: $cmd\n";
			print STDERR "!                                                              !\n";
			print STDERR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
			`mv $logfile $logfile.ERROR`;
			print STDERR "Check $logfile.ERROR for errors...\n\n";
			exit(-1);
		}	

		$no_more_aborts=1;
	}else{
		print STDERR "WARNING: $logfile already exists.\n";
		print STDERR "To rerun in this directory ($output_dir), you need to delete\n";
		print STDERR "the checkpoint/logfile's for the steps you want to repeat.\n";
		$notes="Skipped";
	}

	log_time($date, $time, \@begin_time, \@end_time, "$step_str\_$cmd", $notes);

	$step++;
	return;
}

sub exec_cmd{
	my $exec_string=shift;
	my $log_fname=shift;

	$exec_string=~s/\n//g;
	$exec_string=~s/\t/  /g;
	$exec_string=~s/^\s+//g;

	my ($command_path)=split /\s/, $exec_string;
	my ($command_only, $path)=fileparse($command_path);

	print STDERR "\n";
	print STDERR "***************************************************************\n";
	print STDERR "*  Executing $command_only ...\n";
	print STDERR "***************************************************************\n";

	my $res=`$exec_string 2>&1`;
	
	open(LOG, ">$log_fname.log") || die "Could not open $log_fname.log for writing.\n";
	print LOG "$exec_string\n\n";
	print LOG "$res";
	close(LOG);

	return;
}

sub log_counts{
	my $fname=shift;
	my $comment=shift;
	my $exec_str="$COUNT_NAMES_BIN -n $fname -c $comment -o $output_dir/$COUNTS_LOGNAME";
	`$exec_str`;
}

###############################################################################

# Link fasta file to working directory
my ($name, $path)=fileparse($input_fasta);
print STDERR "Linking $input_fasta -> $output_dir/$name\n";
symlink $input_fasta, "$output_dir/$name";
my $in="$output_dir/$name";
$in=~s/\.fasta$//;

# Link groups file to working directory
my ($group_name, $group_path)=fileparse($groups_file);
print STDERR "Linking $groups_file -> $output_dir/$group_name\n";
symlink $groups_file, "$output_dir/$group_name";
my $group="$output_dir/$group_name";
$group=~s/\.groups//;

# Link 16S alignments and taxa file to working directory
my ($ref16s_name, $ref16s_path)=fileparse($ref_16s_align);
my $reference_name="16S_Reference";
print STDERR "Linking $ref_16s_align -> $output_dir/$reference_name.align\n";
symlink $ref_16s_align, "$output_dir/$reference_name.align";
# Find taxa map file
my $taxa_map=$ref_16s_align;
$taxa_map=~s/\.align$/.tax/;
# Use link as "name"
print STDERR "Linking $taxa_map -> $output_dir/$reference_name.taxa\n";
symlink $taxa_map, "$output_dir/$reference_name.tax";
my $reference_link="$output_dir/$reference_name.align";
my $tax_map="$output_dir/$reference_name.tax";


###############################################################################

execute_mothur_cmd(
	"unique.seqs",
	"fasta=$in.fasta"
);
log_counts("$in.names", "After_1st_Unique");
# Takes
# 	IN.fasta
# Makes 
# 	IN.unique.fasta
#	IN.names

execute_mothur_cmd(
	"align.seqs",
	"candidate=$in.unique.fasta, 
	template=$reference_link, 
	flip=t, 
	processors=$num_proc"
);
# Takes 
# 	IN.unique.fasta
# Makes
# 	IN.unique.align
#	IN.unique.align.report
#	IN.unique.flip.accnos

execute_mothur_cmd(
	"screen.seqs",
	"fasta=$in.unique.align, 
	optimize=start-end, 
	criteria=95,
	name=$in.names,
	group=$group.groups,
	processors=$num_proc"
);
# Takes
# 	GROUP.groups
# Makes
#	IN.unique.good.align
#	IN.unique.bad.accnos
#	IN.good.names
#	GROUP.good.groups
log_counts("$in.good.names", "After_Screening");

execute_mothur_cmd(
	"filter.seqs",
	"fasta=$in.unique.good.align,
	processors=$num_proc"
);
# Makes
#	IN.filter
#	IN.unique.good.filter.fasta


execute_mothur_cmd(
	"unique.seqs",
	"fasta=$in.unique.good.filter.fasta,
	name=$in.good.names"
);
# Makes
#	IN.unique.good.filter.unique.fasta
#	IN.unique.good.filter.names
log_counts("$in.unique.good.filter.names", "After_2nd_Unique");


execute_mothur_cmd(
	"pre.cluster",
	"fasta=$in.unique.good.filter.unique.fasta,
	name=$in.unique.good.filter.names,
	processors=$num_proc"
);
# Makes
#	IN.unique.good.filter.unique.precluster.map
#	IN.unique.good.filter.unique.precluster.fasta
#	IN.unique.good.filter.unique.precluster.names
log_counts("$in.unique.good.filter.unique.precluster.names", "After_Precluster");


execute_mothur_cmd(
	"chimera.uchime",
	"fasta=$in.unique.good.filter.unique.precluster.fasta, 
	name=$in.unique.good.filter.unique.precluster.names, 
	reference=self,
	processors=$num_proc"
);
# Makes
#	IN.unique.good.filter.unique.precluster.denovo.uchime.accnos
#	IN.unique.good.filter.unique.precluster.denovo.uchime.chimeras


execute_mothur_cmd(
	"remove.seqs",
	"accnos=$in.unique.good.filter.unique.precluster.denovo.uchime.accnos, 
	fasta=$in.unique.good.filter.unique.precluster.fasta, 
	name=$in.unique.good.filter.unique.precluster.names,
	group=$group.good.groups"
);
# Makes
# 	IN.unique.good.filter.unique.precluster.pick.names
# 	IN.unique.good.filter.unique.precluster.pick.fasta
#	IN.good.pick.groups
log_counts("$in.unique.good.filter.unique.precluster.pick.names", "After_Chimera_Check");



execute_mothur_cmd(
	"classify.seqs",
	"fasta=$in.unique.good.filter.unique.precluster.pick.fasta, 
	name=$in.unique.good.filter.unique.precluster.pick.names, 
	template=$ref_16s_align,
	taxonomy=$tax_map,
	cutoff=80,
	group=$group.good.pick.groups,
	processors=$num_proc"
);
# Makes
#	IN.unique.good.filter.unique.precluster.pick.REFERENCE.wang.flip.accnos
#	IN.unique.good.filter.unique.precluster.pick.REFERENCE.wang.tax.summary
#	IN.unique.good.filter.unique.precluster.pick.REFERENCE.wang.taxonomy

execute_mothur_cmd(
	"dist.seqs",
	"fasta=$in.unique.good.filter.unique.precluster.pick.fasta, 
	cutoff=$clust_cutoff,
	processors=$num_proc"
);
# Makes
#	IN.unique.good.filter.unique.precluster.pick.dist

execute_mothur_cmd(
	"cluster",
	"column=$in.unique.good.filter.unique.precluster.pick.dist, 
	name=$in.unique.good.filter.unique.precluster.pick.names"
);
# Makes
#	IN.unique.good.filter.unique.precluster.pick.an.sabund
#	IN.unique.good.filter.unique.precluster.pick.an.rabund
#	IN.unique.good.filter.unique.precluster.pick.an.list


execute_mothur_cmd(
	"make.shared",
	"list=$in.unique.good.filter.unique.precluster.pick.an.list, 
	group=$group.good.pick.groups, 
	label=0.03"
);
# Makes
#	IN.unique.good.filter.unique.precluster.pick.an.<per groups>.rabund
#	IN.unique.good.filter.unique.precluster.pick.an.shared

execute_mothur_cmd(
	"classify.otu",
	"taxonomy=$in.unique.good.filter.unique.precluster.pick.$reference_name.wang.taxonomy,
	list=$in.unique.good.filter.unique.precluster.pick.an.list,
	name=$in.unique.good.filter.unique.precluster.pick.names,
	group=$group.good.pick.groups,
	label=0.03"
);
# Makes
# 	IN.unique.good.filter.unique.precluster.pick.an.0.03.cons.taxonomy
#	IN.unique.good.filter.unique.precluster.pick.an.0.03.cons.tax.summary


###############################################################################

# Convert OTU info into Summary Table  
# 	IN.unique.good.filter.unique.precluster.pick.an.shared into Summary Table

my ($date, $time)=format_datetime();
my @sumtab_start_time=time;

my $out_root=$name;
$out_root=~s/\.fasta$//;

my $st_dir="$output_dir/Summary_Tables";
mkdir $st_dir;

my $exec_string="
	$OTU_TO_ST_BIN
		-i $in.unique.good.filter.unique.precluster.pick.an.shared
		-o $st_dir/$out_root.otu
";
exec_cmd($exec_string, "$st_dir/shared_to_summary_table");

# Convert Taxonomy files into Summary Table
#	Will need IN.unique.good.filter.unique.precluster.pick.REFERENCE.wang.taxonomy
#		  IN.unique.good.filter.unique.precluster.pick.names
#		  IN.good.pick.groups

my $exec_string="
	$TAXA_TO_ST_BIN
		-t $in.unique.good.filter.unique.precluster.pick.$reference_name.wang.taxonomy
		-n $in.unique.good.filter.unique.precluster.pick.names
		-g $group.good.pick.groups
		-o $st_dir/$out_root.taxa
";
exec_cmd($exec_string, "$st_dir/taxonomy_to_summary_table");

my @sumtab_end_time=time;

log_time($date, $time, \@sumtab_start_time, \@sumtab_end_time, "generate.summary_tables", "");

###############################################################################

my ($date, $time)=format_datetime();
my @overall_end_time=times;
log_time($date, $time, \@overall_begin_time, \@overall_end_time, "Completion", "Overall");

###############################################################################

print STDERR "done.\n";

