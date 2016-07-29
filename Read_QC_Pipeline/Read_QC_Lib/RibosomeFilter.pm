#!/usr/bin/env perl

###############################################################

package Read_QC_Lib::RibosomeFilter;

=pod

=head1 NAME

Read_QC_Lib::RibosomeFilter

=head1 AUTHOR

Kelvin Li

=head1 DESCRIPTION


=cut

use Carp;
use strict;

###############################################################################

sub new {
	my $this = shift;
	my $class = ref($this) || $this;
	my $self = {};
	bless $self, $class;
	$self->_initialize();
	return $self;
}

###############################################################################

sub _initialize {
	my $self = shift;

	$self->{temporary_directory}="";	# Temporary results
	$self->{input_fasta_path}="";
	$self->{input_fastq_path}="";
	$self->{output_directory}="";
	$self->{processed_fasta}="";
	$self->{executable_path}="";
	$self->{stderr_log}="stderr";
	$self->{stdout_log}="stdout";
	$self->{diagnostic_exec_path}="";
	$self->{remove_list_path}="";

	# Quality Trim specific
	$self->{input_format}="fastq";		# Could be fastq
	$self->{output_format}="fastq";		# Could be fastq
	$self->{database_list}="";
	$self->{alignment_identity_threshold}="";
	$self->{alignment_coverage_threshold}="";
	$self->{alignment_length_threshold}="";
	$self->{STATS_OUT_NAME}="ribosome_filtered.stats";
	$self->{RIBOSOME_FILTERED_OUT_NAME}="ribosome_filtered.fastq";
}

###############################################################################

sub is_null{
	my $val=shift;

	$val=uc($val);	
	$val=~s/\s+//g;
	
	if(($val eq "NA") ||
	   ($val eq "NULL") ||
	   ($val eq "")
	){
		return(1);
	}else{
		return(0);
	}
}

sub set_input_fastq{
	my $self = shift;
	$self->{input_fastq_path}=shift;
}

sub set_temporary_directory{
	my $self = shift;
	$self->{temporary_directory}=shift;
}

sub set_output_directory{
	my $self = shift;
	$self->{output_path}=shift;
}

sub set_executable_path{
	my $self = shift;
	$self->{executable_path}=shift;
}

sub set_alignment_identity_threshold{
	my $self = shift;
	my $val=shift;
	if(is_null($val)){
		$val="";
	};
	$self->{alignment_identity_threshold}=$val;
}

sub set_alignment_coverage_threshold{
	my $self = shift;
	my $val=shift;
	if(is_null($val)){
		$val="";
	};
	$self->{alignment_coverage_threshold}=$val;
}

sub set_alignment_length_threshold{
	my $self = shift;
	my $val=shift;
	if(is_null($val)){
		$val="";
	};
	$self->{alignment_length_threshold}=$val;
}

sub set_database_list{
	my $self = shift;
	$self->{database_list}=shift;
}

sub get_processed_fastq{
	my $self = shift;
	return($self->{processed_fastq});
}

sub get_output_format{
        my $self = shift;
        return($self->{output_format});
}

sub get_input_format{
        my $self = shift;
        return($self->{input_format});
}

###############################################################################

sub check_variable{
	my $var=shift;
	my $varname=shift;
	if($var eq ""){
		print STDERR "WARNING: $varname is undefined.\n";
		return(1);
	}else{
		return(0);
	}
}

###############################################################################

sub execute_analysis{
	my $self=shift;
	
	my $exec_path=$self->{executable_path};
	my $output_path=$self->{output_path};
	my $fastq_path=$self->{input_fastq_path};	

	# Confirm necessary variables are set
	my $err=0;
	$err+=check_variable($exec_path, "Execution Path");
	$err+=check_variable($output_path, "Result Path");
	$err+=check_variable($fastq_path, "Input FASTQ Path");
	die "Variables undefined." unless !$err;

	# Make output directory if necessary
	if(!(-e $output_path)){
		mkdir $output_path;
	}

	# Construct option string for thresholds
	my $thresholds="";

	my $idn_thres=$self->{alignment_identity_threshold};
	my $cov_thres=$self->{alignment_coverage_threshold};
	my $len_thres=$self->{alignment_length_threshold};

	if(!check_variable($idn_thres, "Alignment Identity Threshold")){
		$thresholds.="-i $idn_thres ";
	}
	if(!check_variable($idn_thres, "Alignment Coverage Threshold")){
		$thresholds.="-c $cov_thres ";
	}
	if(!check_variable($len_thres, "Alignment Length Threshold")){
		$thresholds.="-l $len_thres ";
	}

	# Get database list
	my $database_list=$self->{database_list};
	die "Database not specified." unless !check_variable($database_list, "Database List");

	# Construct execution string
	my $execute_analysis_string=
		"$exec_path " .
		"-f $fastq_path " .
		"-dbs $database_list " .
		"$thresholds " .
		"-keep_tmp_files " .
		"-out_dir $output_path";

	# Execute
	print STDERR "Executing: $execute_analysis_string\n";
	#my $res=`$execute_analysis_string`;
	system($execute_analysis_string);

	# Keep track of where final output fastq file is
	my @output_dir_list=split "\n", `find $output_path | grep _nonrrna.fq`;
	if($#output_dir_list!=0){
		die "Error, could not find exactly one _nonrrna.fq file in output directory\n";
	}
	$self->{processed_fastq}=$output_dir_list[0];
	
	return;
}	

sub execute_diagnostics{
	my $self=shift;
	my $diag_exec=$self->{diagnostic_exec_path};
	if($diag_exec eq ""){
		print STDERR "No diagnostics specified.\n";
	}else{
		# Insert diagnostics here
	}
	return;
}

sub perform_qc{

	my $self=shift;

	# 1. Run analysis
	print STDERR "\nRunning Analysis:\n";
	$self->execute_analysis();

	print STDERR "\nDone.\n";

	return;
}

###############################################################################

sub DESTROY {
    my $self = shift;
}

###############################################################################

1;

