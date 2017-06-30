#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../Annotation_Lib";
use Getopt::Std;
use FileHandle;
use Config::IniFiles;
use POSIX; # for ceil function
use Cwd;

use vars qw ($opt_l $opt_p $opt_o $opt_r $opt_s $opt_e);

getopts("l:p:o:r:s:");

my $usage = "
usage:
    $0

	-l <Filename for List of Query FASTA files>
	
	-p <parameter .ini file>
	-o <output directory>
	
	[-r <input fasta root directory>]
	[-s <Offset, default=1>,<Multiplier, default=1>]
	[-e (execute/launch, or just produce commands in shell script)]

	This script will load the fasta list, where each fasta file represents a single
	sample/library, and run the specified aligment tool.  For each sample/library
	a directory in the output directory will be created and the alignment output
	will be saved there, and converted into a composite alignment output file.

	If the -s option is used, then only a subset of the FASTA files 
        will be processed. If offset=1, and multipler=1, then
        items 1,2,3, ... N sample will be processed, i.e. all  If offset=3 and
        mulitplier=2, then 3,5,7, ... N sample will be processed.
        Offset can be thought of the job ID starting from 1, and multiplier
        is the number of processes to run in parallel. 

	INPUT:

	The -l option specifies the sample ID and FASTA file in the following format:

		<sample_id> \\t <fasta file (relative path)> \\n

	OUTPUT:

	In <output directory>/<sample id>:
	
	1.) Alignment tool output (gzip'd)

	2.) Composite output:
	
		<query id> <subject id> <percent composite> 



";

###############################################################################

if(
	!defined($opt_l) ||
	!defined($opt_p) ||
	!defined($opt_o)
){
	die $usage;
}

my $FASTA_List_Fn=$opt_l;
my $Parameter_Fn=$opt_p;
my $Output_Directory=$opt_o;
my $Input_FASTA_Root_Directory=$opt_r;
my $Offset;
my $Multiplier;
my $Launch;

if(defined($opt_s)){
	($Offset, $Multiplier)=split ",", $opt_s;
}else{
	($Offset, $Multiplier)=(1,1);
}

if(defined($opt_e)){
	$Launch=1;
}else{
	$Launch=0;
}

print STDERR "FASTA List: $FASTA_List_Fn\n";
print STDERR "Parameter Filename: $Parameter_Fn\n";
print STDERR "Output Directory: $Output_Directory\n";
print STDERR "\n";
print STDERR "Input FASTA Root Directory: $Input_FASTA_Root_Directory\n";
print STDERR "Offset: $Offset  Multiplier: $Multiplier\n";
print STDERR "Launch: $Launch\n";
print STDERR "\n";

###############################################################################

sub load_file_list{
        my $list=shift;

        system("dos2unix $list");

        open(IN_FH, "<$list") || die "Could not opne $list\n";
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

sub get_config_val{
	my $config=shift;
	my $group=shift;
	my $field=shift;

	my $value=$config->val($group, $field);
	
	# Dereference variable pointing elsewhere
	if($value=~/\[(\w+)\]:(\w+)/){
		$value=get_config_val($config, $1, $2);
		print STDERR "Resolved [$1]:$2 to $value\n";
	}
		
	return($value);

}

###############################################################################

sub log_timing_stats{
        my $time_rec_begin_ref=shift;
        my $time_rec_end_ref=shift;
        my $time_output_fn=shift;
        my $name=shift;

        my $start=0;
        my $end=0;

        for(my $i=0; $i<4; $i++){
                $start+=${$time_rec_begin_ref}[$i];
                $end+=${$time_rec_end_ref}[$i];
        }

        my $tot_time=$end-$start;

        open(FH, ">>$time_output_fn") || die "Could not open $time_output_fn for appending.\n";

        print FH "$name\t$tot_time\n";

        close(FH);
}

###############################################################################

sub align{
	my $query_fasta_file=shift;
	my $output_directory=shift;
	my $launch=shift;
	my $config=shift;


	my $cfg_bl_program=get_config_val($config, "blast_options", "program");
	my $cfg_bl_eval=get_config_val($config, "blast_options", "eval");
	my $cfg_bl_dbsize=get_config_val($config, "blast_options", "dbsize");
	my $cfg_bl_num_threads=get_config_val($config, "blast_options", "num_threads");
	my $cfg_bl_tool=get_config_val($config, "blast_options", "aligner");

	my $aln_cmd;
	my $cfg_bl_db;

	my $cwd=cwd();
	my $out_dir_w_wd=$output_directory;
	if($output_directory=~/^\.\//){
		my $clean_od=$output_directory;
		$clean_od=~s/\.\///;
		$out_dir_w_wd="\$working_directory/$clean_od";
	}

	if($cfg_bl_tool eq "ncbi_blast"){

		my $cfg_bl_program_bin=get_config_val($config, "blast_options", "program");
		$cfg_bl_db=get_config_val($config, "blast_options", "blast_db");

		$aln_cmd="
		$cfg_bl_program_bin
		-query $query_fasta_file
		-db $cfg_bl_db
		-out $output_directory/blastx.out
		-evalue $cfg_bl_eval
		-outfmt 6
		-dbsize $cfg_bl_dbsize
		-num_threads $cfg_bl_num_threads
		";

	}elsif($cfg_bl_tool eq "diamond"){

		my $cfg_bl_diamond_bin=get_config_val($config, "blast_options", "diamond_bin");
		my $cfg_bl_program=get_config_val($config, "blast_options", "program");
		$cfg_bl_db=get_config_val($config, "blast_options", "diamond_db");

		$aln_cmd="
		$cfg_bl_diamond_bin
		$cfg_bl_program
		--query $query_fasta_file
		--db $cfg_bl_db
		--out $out_dir_w_wd/blastx.out
		--evalue $cfg_bl_eval
		--outfmt 6
		--dbsize $cfg_bl_dbsize
		--threads $cfg_bl_num_threads
		";

	}

	$aln_cmd=~s/\s+/ /g;
	print STDERR "Command:\n";
	print STDERR "$aln_cmd\n";

	# Compuate percent composite identity
	my $cfg_perc_comp_id_bin=get_config_val($config, "tools", "perc_comp_id_bin");
	my $cmpos_cmd="
		perl $cfg_perc_comp_id_bin
		-q $query_fasta_file 
		-p $out_dir_w_wd/blastx.out 
	";
	$cmpos_cmd=~s/\s+/ /g;
        print STDERR "Command:\n";
        print STDERR "$cmpos_cmd\n";

	# Output shell script
	my $algn_cmd_shsc="$output_directory/align.csh";
	open(FH, ">$algn_cmd_shsc") || die "Could not open $algn_cmd_shsc\n";
	print FH "#!/bin/csh\n";
	print FH "\n";
	print FH "set working_directory=$cwd\n";
	print FH "\n";
	print FH "echo 'Starting alignment of:'\n";
	print FH "echo '\t' $query_fasta_file\n";
	print FH "echo 'against:'\n";
	print FH "echo '\t' $cfg_bl_db...\n";
	print FH "\n";
	print FH "$aln_cmd\n";
	print FH "\n";
	print FH "echo 'alignment finished.'\n";
	print FH "\n";
	print FH "echo 'Computing Percent Composite Identity (PCI).'\n";
	print FH "\n";
	print FH "$cmpos_cmd\n";
	print FH "\n";
	print FH "echo 'PCI compute completed.'\n";
	print FH "\n";
	print FH "echo 'done.'\n";
	close(FH);

	`chmod +x $algn_cmd_shsc`;

	# Execute
	if($launch){
		system("$algn_cmd_shsc");
	}
	
}

###############################################################################
###############################################################################

my $cfg=Config::IniFiles->new(-file => $Parameter_Fn);

my $ini_global_scratch_dir=get_config_val($cfg, "GLOBAL", "scratch_dir");
print STDERR "Scratch directory: $ini_global_scratch_dir\n";

my $file_arr_ref=load_file_list($FASTA_List_Fn);
my $num_records=$#{$file_arr_ref}+1;

print STDERR "Num records in FASTA list: $num_records\n";

# Try to make output dir
make_dir($Output_Directory);

my $offset_str_width=ceil(log($Multiplier)/log(10));
my $offset_str=sprintf("%0$offset_str_width"."i", $Offset);
my $timing_log_fn="$Output_Directory/timing_log.$offset_str.tsv";

$Offset--; # Start offset from 1 less than specified, so we start from 0, instead of 1
for(my $idx=$Offset; $idx<$num_records; $idx+=$Multiplier){

	print STDERR "\nWorking on record: $idx.\n";

	my @rec=split "\t", ${$file_arr_ref}[$idx];
	my ($sample_id, $fname)=@rec;

	print STDERR "Target Sample ID: $sample_id\n";
	print STDERR "Target Filename: $fname\n";

	my $sample_output_dir="$Output_Directory/$sample_id";
	make_dir($sample_output_dir);

	if($Input_FASTA_Root_Directory ne ""){
		$fname="$Input_FASTA_Root_Directory/$fname";
	}

	align($fname, $sample_output_dir, $Launch, $cfg);

}

###############################################################################

print STDERR "Done.\n";

