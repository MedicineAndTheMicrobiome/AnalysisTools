#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Std;
use vars qw($opt_l $opt_p $opt_s $opt_r $opt_c $opt_o);
use File::Basename;
use POSIX; # for ceil function

use Read_QC_Lib::Dust;
use Read_QC_Lib::AdapterTrim;
use Read_QC_Lib::ContaminantFilter;
use Read_QC_Lib::QualityTrim;
use Read_QC_Lib::RibosomeFilter;
use Read_QC_Lib::Subsample;

use Read_QC_Lib::Config::IniFiles;

my $command_list = "contam,dust,qvtrim,adapt,ribo";

getopts("l:p:s:r:c:o:");
my $usage = "usage: 

$0 

	-l <List of reads to process>
	-p <Parameter file>
	-o <Output directory>
	[-r <input fastq root directory>]
	[-s <Offset, default=1>,<Multiplier, default=1>]
	[-c <Command list, default=\"$command_list\"]
	
	This script will run the list of read processing steps
	specified in the order of the command list.

	The format of the list is:

	<library name> <for fastq path> <rev fastq path> <for adapter> <rev adapter>

	Each column is tab separated.
	
	If the -s option is used, then only a subset of the reads
	will be processed. If offset=1, and multipler=1, then
	items 1,2,3, ... N sample will be processed.  If offset=3 and
	mulitplier=2, then 3,5,7, ... N sample will be processed.
	Offset can be thought of the job ID starting from 1, and multiplier
	is the number of processes to run in parallel.  

	The output directory (-o), is the directory where each of the processed
	libraries will be saved.

	Possible Commands:

		subsample:	Subsample
		contam:		(Host) Contaminant Screen
		dust:		Low Complexity Filter
		qvtrim:		Quality Value
		barcode_trim:	Barcode Trimming
		seq_adapt_trim:	Adapter Trimming
		primer_trim:	Primer Trimming
		ribo:		Ribosomal Screening

	An example of a quick QC run would be:
		subsample,contam,dust,qvtrim,adapt,ribo

	If the -r option is specified, then that path will be pre-pended to the forward and
	reverse fastq paths if they are relative.  i.e. they do not begin with a /.
	
	
";

if(!defined($opt_l)||!defined($opt_p)||!defined($opt_o)){
	die $usage;
}

my $read_list=$opt_l;
my $parameter_fname=$opt_p;
my $output_dir=$opt_o;
my $subparam="1,1";
my $input_fastq_root="";


if(defined($opt_s)){
	$subparam=$opt_s;	
}

if(defined($opt_c)){
	$command_list=$opt_c;
}

if(defined($opt_r)){
	$input_fastq_root=$opt_r;
}

my ($offset, $multiplier)=split ",", $subparam;

my @command_arr=split ",", $command_list;

###############################################################################

print STDERR "\n";
print STDERR "Read List: $read_list\n";
print STDERR "Parameters: $parameter_fname\n";
print STDERR "Subset Paramters: Offset=$offset Multiplier=$multiplier\n";
print STDERR "Input fastq root director: $input_fastq_root\n";

print STDERR "Command List:\n";
foreach my $command(@command_arr){
	print STDERR "\t$command\n";
}

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
		print STDERR "$output_dir exists.\n";
	}
	if(!(-e $dir)){
		die "Could not make or find $dir\n";
	}
}

###############################################################################

sub convert_fastq_to_fasta{
	my $fastq_path=shift;
	my $config=shift;

	# Generate output fasta file name
	my $fasta_path=$fastq_path;
	$fasta_path=~s/\.fastq.gz$/.fasta/;
	$fasta_path=~s/\.fastq$/.fasta/;
	$fasta_path=~s/\.fq.gz$/.fasta/;
	$fasta_path=~s/\.fq$/.fasta/;

	my $path=$config->val("FASTQ_Tools", "Path");
	my $bin=$config->val("FASTQ_Tools", "SplitFASTQ_bin");

	print STDERR "Running Convert FASTQ to FASTA: $bin\n";
	my $res=`$path/$bin -f $fastq_path -a $fasta_path`;

	return($fasta_path);
}

sub filter_fastq_by_list{
	my $fastq_path=shift;
	my $remove_list_path=shift;
	my $dest_path=shift;
	my $config=shift;	

	# Generate output fastq file name
	my $out_fastq=$fastq_path;
	$out_fastq=~s/\.fastq$//;
	$out_fastq=~s/\.fq$//;
	$out_fastq.=".filt.fastq";

	my ($name, $path)=fileparse($out_fastq);
	$out_fastq="$dest_path/$name";

	my $path=$config->val("FASTQ_Tools", "Path");
	my $bin=$config->val("FASTQ_Tools", "Exclude_by_ID_bin");

	print STDERR "Running Filter FASTQ by List: $bin\n";
	my $res=`$path/$bin -f $fastq_path -l $remove_list_path -o $out_fastq`;

	return($out_fastq);
}

sub merge_paired_fastq{
	my $forward_path=shift;
	my $reverse_path=shift;
	my $output_path=shift;
	my $config=shift;

	my $path=$config->val("FASTQ_Tools", "Path");
	my $bin=$config->val("FASTQ_Tools", "Merge_Paired_FASTQ_bin");
	my $tmp_dir=$config->val("GLOBAL", "scratch_dir");
	
	print STDERR "Running Merge Paired FASTQ: $bin\n";
	my $res=`$path/$bin -f $forward_path -r $reverse_path -o $output_path -t $tmp_dir`;

	return($output_path);
}

sub log_fasta_stats{
	my $target_fastq_fn=shift;
	my $stat_output_fn=shift;
	my $name=shift;
	my $config=shift;

	my $path=$config->val("FASTQ_Tools", "Path");
	my $bin=$config->val("FASTQ_Tools", "Report_FASTQ_Stats_bin");
	
	print STDERR "Running Report FASTQ Stats: $bin\n";

	if(-e $stat_output_fn){
		my $res=`$path/$bin -f $target_fastq_fn -h -a \"$name\" >> $stat_output_fn`;
	}else{
		my $res=`$path/$bin -f $target_fastq_fn -a \"$name\" > $stat_output_fn`;
	}

	return;
}

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


sub execute_module{
	my $command=shift;	
	my $working_directory=shift;
	my $input_fastq_path=shift;
	my $output_directory=shift;
	my $config=shift;
	my $run_specific_param_hash_ref=shift;
	
	my $VERBOSE=1;
	if($VERBOSE){
		print STDERR "Command: $command\n";
		print STDERR "Working Directory: $working_directory\n";
		print STDERR "Input FASTQ: $input_fastq_path\n";
		print STDERR "Output Dir: $output_directory\n";
	}

	my $global_scratch=$config->val("GLOBAL", "scratch_dir");
	my $mod;

	# Set up the QC module based on the command
	if($command eq "contam"){

		$mod=new Read_QC_Lib::ContaminantFilter;
		$mod->set_executable_path($config->val("contaminant_screen", "bin"));
		$mod->set_reference_srprism($config->val("contaminant_screen", "reference_srprism"));
		$mod->set_reference_bitmask($config->val("contaminant_screen", "reference_bitmask"));
		$mod->set_reference_blastdb($config->val("contaminant_screen", "reference_blastdb"));

	}elsif($command eq "dust"){

		$mod=new Read_QC_Lib::Dust;
		$mod->set_executable_path($config->val("dust_screen", "bin"));
		$mod->set_cutoff($config->val("dust_screen", "cutoff"));

	}elsif($command eq "qvtrim"){

		$mod=new Read_QC_Lib::QualityTrim;
		$mod->set_trim_executable_path($config->val("quality_trim", "bin_trim"));
		$mod->set_filt_executable_path($config->val("quality_trim", "bin_filt"));

		$mod->set_trim_quality_threshold($config->val("quality_trim", "trim_quality_threshold"));
		$mod->set_trim_minimum_length($config->val("quality_trim", "trim_minimum_length"));
		$mod->set_filt_quality_threshold($config->val("quality_trim", "filt_quality_threshold"));
		$mod->set_filt_percent_above($config->val("quality_trim", "filt_percent_above"));

	}elsif($command eq "barcode_trim"){
		
		my $adapt_seq_F=${$run_specific_param_hash_ref}{"barcode.forward"};
		my $adapt_seq_R=${$run_specific_param_hash_ref}{"barcode.reverse"};
		if(!defined($adapt_seq_F) || !defined($adapt_seq_R)){
			die "Barcode sequence not specified in run specific parameter hash: barcode.forward / barcode.reverse\n";
		}

		$mod=new Read_QC_Lib::AdapterTrim;
		$mod->set_executable_path($config->val("adapter_trim", "bin"));
		$mod->add_adapter_sequence($adapt_seq_F);
		$mod->add_adapter_sequence($adapt_seq_R);

	}elsif($command eq "seq_adapt_trim"){
		
		my $seq_adapt_F=$config->val("seq_adapt_trim", "forward");
		my $seq_adapt_R=$config->val("seq_adapt_trim", "reverse");

		if(!defined($seq_adapt_F) || !defined($seq_adapt_R)){
			die "Sequencing adapter sequence not specified in configuration file.\n";
		}

		$mod=new Read_QC_Lib::AdapterTrim;
		$mod->set_executable_path($config->val("adapter_trim", "bin"));
		$mod->add_adapter_sequence($seq_adapt_F);
		$mod->add_adapter_sequence($seq_adapt_R);

	}elsif($command eq "primer_trim"){
		
		my $primer_seq_F=$config->val("primer_trim", "forward");
		my $primer_seq_R=$config->val("primer_trim", "reverse");

		if(!defined($primer_seq_F) || !defined($primer_seq_R)){
			die "Primer sequence not specified in configuration file.\n";
		}

		$mod=new Read_QC_Lib::AdapterTrim;
		$mod->set_executable_path($config->val("adapter_trim", "bin"));
		$mod->add_adapter_sequence($primer_seq_F);
		$mod->add_adapter_sequence($primer_seq_R);

	}elsif($command eq "ribo"){

		$mod=new Read_QC_Lib::RibosomeFilter;	
		$mod->set_executable_path($config->val("ribosomal_screen", "bin"));
		$mod->set_alignment_identity_threshold($config->val("ribosomal_screen", "alignment_identity_threshold"));
		$mod->set_alignment_coverage_threshold($config->val("ribosomal_screen", "alignment_coverage_threshold"));
		$mod->set_alignment_length_threshold($config->val("ribosomal_screen", "alignment_length_threshold"));
		$mod->set_database_list($config->val("ribosomal_screen", "database_list"));

	}elsif($command eq "subsample"){

		$mod=new Read_QC_Lib::Subsample;
		$mod->set_executable_path($config->val("subsample", "bin"));
		$mod->set_sample_size($config->val("subsample", "sample_size"));

	}else{
		die "Undefined command: $command\n";
	}

	$mod->set_temporary_directory($global_scratch);
	$mod->set_output_directory($output_directory);

	# Make sure we have the right input format
	my $input_format=$mod->get_input_format();
	if($input_format eq "fastq"){
		if($VERBOSE){print STDERR "Setting fastq to $input_fastq_path\n";}
		$mod->set_input_fastq($input_fastq_path);
	}elsif($input_format eq "fasta"){
		my $fasta_path=convert_fastq_to_fasta($input_fastq_path, $config);
		if($VERBOSE){print STDERR "Setting fasta to $fasta_path\n";}
		$mod->set_input_fasta($fasta_path);
	}

	# Run the program
	$mod->perform_qc();

	# Make sure we generate a fastq file
	my $fastq_location;
	my $result_location;
	my $output_format=$mod->get_output_format();
	if($output_format eq "fastq"){
		$result_location=$mod->get_processed_fastq();
	}elsif($output_format eq "remove list"){
		my $remove_list_path=$mod->get_remove_list_path();
		my ($name, $path)=fileparse($remove_list_path);
		$result_location=filter_fastq_by_list($input_fastq_path, $remove_list_path, $path, $config);
	}

	return($result_location);
}



###############################################################################

my $cfg=Config::IniFiles->new(-file => $parameter_fname);

my $ini_global_scratch_dir=$cfg->val("GLOBAL", "scratch_dir");
print STDERR "Scratch directory: $ini_global_scratch_dir\n";

my $file_hash_ref=load_file_list($read_list);
my $num_commands=$#command_arr+1;

###############################################################################

# Try to make output dir
make_dir($output_dir);

my $offset_str_width=ceil(log($multiplier)/log(10));
my $offset_str=sprintf("%0$offset_str_width"."i", $offset);

my $timing_log_fn="$output_dir/timing_log.$offset_str.tsv";
my $fastq_stat_log_fn="$output_dir/fastq_qc_log.$offset_str.tsv";

my $UNSPECIFIED="[Unspecified]";
my $NA="[Not Applicable]";

my $num_records=$#{$file_hash_ref}+1;
$offset--; # Start offset from 1 less than specified, so we start from 0, instead of 1
for(my $idx=$offset; $idx<$num_records; $idx+=$multiplier){

	print STDERR "\nWorking on record: $idx.\n";
	my ($lib_name, $for_path, $rev_path, $for_adapt, $rev_adapt)= split "\t", ${$file_hash_ref}[$idx];

	if($input_fastq_root ne ""){
		if(substr($for_path, 0 , 1) ne "/"){
			$for_path="$input_fastq_root/$for_path";
		}
		if(substr($rev_path, 0 , 1) ne "/"){
			$rev_path="$input_fastq_root/$rev_path";
		}
	}

	if($rev_path eq ""){
		$rev_path = $UNSPECIFIED;
		$rev_adapt = $NA;
	}

	my @paired_time_rec_begin=times;

	print STDERR "\n";
	print STDERR "###############################################################################\n";
	print STDERR "# Library Name: $lib_name\n";
	print STDERR "###############################################################################\n";
	print STDERR "\n";
	print STDERR "  Forward Path: $for_path\n";
	print STDERR "  Forward Adapter: $for_adapt\n";
	print STDERR "  Reverse Path: $rev_path\n";
	print STDERR "  Reverse Adapter: $rev_adapt\n";
	print STDERR "\n";

	
	# Make directory for library
	my $library_dir="$output_dir/$lib_name";
	if(-e "$library_dir/completed."){
		print STDERR "It appears $library_dir has already successfully completed.  Skipping.\n";
		next;
	}
	make_dir($library_dir);

	# Set up has for library/iteration specific parameters
	my %run_spec_params;
	$run_spec_params{"barcode.forward"}=$for_adapt;
	$run_spec_params{"barcode.reverse"}=$rev_adapt;

	# Set up loop to do both directions individually
	my %input_reads_path_hash;
	my %output_reads_path_hash;

	my @dir_list=("F", "R");
	if($rev_path eq $UNSPECIFIED){
		@dir_list=("F");
	}

	$input_reads_path_hash{"F"}=$for_path;
	$input_reads_path_hash{"R"}=$rev_path;
	$output_reads_path_hash{"F"}=$for_path;
	$output_reads_path_hash{"R"}=$rev_path;

	# Do each direction independently
	foreach my $direction (@dir_list){

		make_dir("$library_dir/$direction");

		my ($name, $path)=fileparse($input_reads_path_hash{$direction});

		my $suc=symlink("$path/$name", "$library_dir/$direction/$name");
		my $current_fastq_fn="$library_dir/$direction/$name";

		log_fasta_stats($current_fastq_fn, $fastq_stat_log_fn,
			"$lib_name\t$direction\traw", $cfg);

		# Step through all the QC steps
		for(my $cmd_idx=0; $cmd_idx<$num_commands; $cmd_idx++){

			my $cur_command=$command_arr[$cmd_idx];
			my $cmd_tag="$cmd_idx\_$cur_command";

			print STDERR "\n";
			print STDERR "***************************************\n";
			print STDERR "* Performing command: $cmd_tag\n";
			print STDERR "***************************************\n";
			print STDERR "\n\n";

			make_dir("$library_dir/$direction/$cmd_tag");

			my @time_rec_begin=times;

			$output_reads_path_hash{$direction}=
				execute_module(
					$cur_command,
					"$ini_global_scratch_dir/$cmd_tag",
					$current_fastq_fn,
					"$library_dir/$direction/$cmd_tag",
					$cfg,
					\%run_spec_params
				);

			my @time_rec_end=times;

			# Set output from current as input of next step
			$current_fastq_fn=$output_reads_path_hash{$direction};


			# Log the fastq stats
			log_fasta_stats($current_fastq_fn, $fastq_stat_log_fn,
				"$lib_name\t$direction\t$cur_command", $cfg);
			log_timing_stats(\@time_rec_begin, \@time_rec_end, $timing_log_fn, 
				"$lib_name\t$direction\t$cur_command");

		}	

	}


	if($rev_path ne $UNSPECIFIED){
		# Join both directions after all processing has finished
		merge_paired_fastq(
			$output_reads_path_hash{"F"},
			$output_reads_path_hash{"R"},
			"$library_dir/$lib_name",
			$cfg
		);

		log_fasta_stats("$library_dir/$lib_name.for.frag.fastq", $fastq_stat_log_fn,
			"$lib_name\tfor_frag\tmerged", $cfg);
		log_fasta_stats("$library_dir/$lib_name.rev.frag.fastq", $fastq_stat_log_fn,
			"$lib_name\trev_frag\tmerged", $cfg);
		log_fasta_stats("$library_dir/$lib_name.paired.fastq", $fastq_stat_log_fn,
			"$lib_name\tpaired\tmerged", $cfg);

		my @paired_time_rec_end=times;
		log_timing_stats(\@paired_time_rec_begin, \@paired_time_rec_end, $timing_log_fn, 
			"$lib_name\tPaired\tcombined");


	}else{
		# The Forward file is the final file
		# Make it a relative
		my $src_path=$output_reads_path_hash{"F"};
		$src_path=~s/^$library_dir\///;
		symlink $src_path, "$library_dir/$lib_name.qcd.fasta";

		log_fasta_stats("$library_dir/$lib_name.qcd.fasta", $fastq_stat_log_fn,
			"$lib_name\tForward\tfinal", $cfg);

		my @paired_time_rec_end=times;
		log_timing_stats(\@paired_time_rec_begin, \@paired_time_rec_end, $timing_log_fn, 
			"$lib_name\tIndividual\tfinal");

	}

	# Make a file when all is done.
	`touch $library_dir/completed.`;


}

my @begin=(0,0,0,0);
my @time_final_end=times;
my ($name, $path)=fileparse($read_list);
log_timing_stats(\@begin, \@time_final_end, $timing_log_fn, "Total Execution\t$name\t");

###############################################################################

print STDERR "Done.\n";

