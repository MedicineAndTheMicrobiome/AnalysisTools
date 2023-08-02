#!/usr/bin/env perl

###############################################################

package Read_QC_Lib::AbsoluteTrim;

=pod

=head1 NAME

Read_QC_Lib::AbsoluteTrim

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
	$self->{keep_start}="1";
	$self->{keep_end}="Inf";
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

sub set_executable_path{
	my $self = shift;
	$self->{executable_path}=shift;
}

sub set_trim_start_keep{
	my $self = shift;
	$self->{keep_start}=shift;
}

sub set_trim_end_keep{
	my $self = shift;
	$self->{keep_end}=shift;
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
	
	my $trim_exec_path=$self->{executable_path};
	my $output_path=$self->{output_path};
	my $fastq_path=$self->{input_fastq_path};	

	my $trim_keep_start=$self->{keep_start};
	my $trim_keep_end=$self->{keep_end};

	# Confirm necessary variables are set
	my $err=0;
	$err+=check_variable($trim_exec_path, "Trimmer Execution Path");
	$err+=check_variable($output_path, "Result Path");
	$err+=check_variable($fastq_path, "Input FASTQ Path");
	$err+=check_variable($trim_keep_start, "Trimmer Start Keep");
	$err+=check_variable($trim_keep_end, "Trimmer End Keep");
	die "Variables undefined." unless !$err;

	my $end_opt_str;
	if(uc($trim_keep_end) eq "INF"){
		$end_opt_str="";
	}else{
		$end_opt_str="-l $trim_keep_end ";
	}

	# Make output directory if necessary
	if(!(-e $output_path)){
		mkdir $output_path;
	}

	# Construct execution string
	my $execute_analysis_string=
		"cat $fastq_path | $trim_exec_path " .
		"  -f $trim_keep_start " .
		$end_opt_str .
		"> $output_path/abs_trimmed.fastq";

	# Execute
	print STDERR "Executing: $execute_analysis_string\n";
	my $res=`$execute_analysis_string`;

	$self->{processed_fastq}="$output_path/abs_trimmed.fastq";
	
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

