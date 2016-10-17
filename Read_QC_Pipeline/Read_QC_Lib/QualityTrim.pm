#!/usr/bin/env perl

###############################################################

package Read_QC_Lib::QualityTrim;

=pod

=head1 NAME

Read_QC_Lib::QualityTrim

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
	$self->{quality_threshold}="";
	$self->{minimum_length}="";
	$self->{STATS_OUT_NAME}="qual.out.stats";
	$self->{QVTRIM_OUT_NAME}="qv_trimmed.fastq";
}

###############################################################################

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

sub set_trim_executable_path{
	my $self = shift;
	$self->{trim_executable_path}=shift;
}

sub set_filt_executable_path{
	my $self = shift;
	$self->{filt_executable_path}=shift;
}

sub set_trim_quality_threshold{
	my $self = shift;
	$self->{trim_quality_threshold}=shift;
}

sub set_trim_minimum_length{
	my $self = shift;
	$self->{trim_minimum_length}=shift;
}

sub set_filt_quality_threshold{
	my $self = shift;
	$self->{filt_quality_threshold}=shift;
}

sub set_filt_percent_above{
	my $self = shift;
	$self->{filt_percent_above}=shift;
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
	
	my $trim_exec_path=$self->{trim_executable_path};
	my $filt_exec_path=$self->{filt_executable_path};
	my $output_path=$self->{output_path};
	my $fastq_path=$self->{input_fastq_path};	

	my $trim_quality_threshold=$self->{trim_quality_threshold};
	my $trim_minimum_length=$self->{trim_minimum_length};
	my $filt_quality_threshold=$self->{filt_quality_threshold};
	my $filt_percent_above=$self->{filt_percent_above};

	# Confirm necessary variables are set
	my $err=0;
	$err+=check_variable($trim_exec_path, "Trimmer Execution Path");
	$err+=check_variable($filt_exec_path, "Filter Execution Path");
	$err+=check_variable($output_path, "Result Path");
	$err+=check_variable($fastq_path, "Input FASTQ Path");
	$err+=check_variable($trim_quality_threshold, "Trimmer Quality Threshold");
	$err+=check_variable($trim_minimum_length, "Trimer Minimum Length");
	$err+=check_variable($filt_quality_threshold, "Filter Quality Threshold");
	$err+=check_variable($filt_percent_above, "Filter Percent Above");
	die "Variables undefined." unless !$err;

	# Make output directory if necessary
	if(!(-e $output_path)){
		mkdir $output_path;
	}

	# Construct execution string
	my $execute_analysis_string=
		"$trim_exec_path " .
		"  -i $fastq_path " .
		"  -t $trim_quality_threshold -l $trim_minimum_length -v | " .
		"$filt_exec_path " .
		"  -q $filt_quality_threshold -p $filt_percent_above -v " .
		"> $output_path/$self->{QVTRIM_OUT_NAME}";

	# Execute
	print STDERR "Executing: $execute_analysis_string\n";
	my $res=`$execute_analysis_string`;

	$self->{processed_fastq}="$output_path/$self->{QVTRIM_OUT_NAME}";
	
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

