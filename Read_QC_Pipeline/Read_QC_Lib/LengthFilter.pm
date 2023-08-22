#!/usr/bin/env perl

###############################################################

package Read_QC_Lib::LengthFilter;

=pod

=head1 NAME

Read_QC_Lib::LengthFilter

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
	$self->{input_format}="fasta";		# Could be fastq
	$self->{output_format}="remove list";	# Could be fastq
	$self->{keep_start}="1";
	$self->{keep_end}="Inf";
}

###############################################################################

sub set_input_fasta{
	my $self = shift;
	$self->{input_fasta_path}=shift;
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

sub set_min_length{
	my $self = shift;
	$self->{min_length}=shift;
}

sub set_max_length{
	my $self = shift;
	$self->{max_length}=shift;
}

sub get_remove_list_path{
	my $self = shift;
	return($self->{remove_list_path});
}

sub get_processed_fasta{
	my $self = shift;
	return($self->{processed_fasta});
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
	my $fasta_path=$self->{input_fasta_path};	

	my $min_length=$self->{min_length};
	my $max_length=$self->{max_length};

	# Confirm necessary variables are set
	my $err=0;
	$err+=check_variable($exec_path, "Length Filter Execution Path");
	$err+=check_variable($output_path, "Result Path");
	$err+=check_variable($fasta_path, "Input FASTA Path");
	$err+=check_variable($min_length, "Min Length");
	$err+=check_variable($max_length, "Max Length");
	die "Variables undefined." unless !$err;

	# Make output directory if necessary
	if(!(-e $output_path)){
		mkdir $output_path;
	}

	# Construct execution string
	my $execute_analysis_string=
		"$exec_path " .
		"  -f $fasta_path " .
		"  -n $min_length " .
		"  -x $max_length " .
		"  -r " .
		"> $output_path/length.exclusion.lst";

	# Execute
	print STDERR "Executing: $execute_analysis_string\n";
	my $res=`$execute_analysis_string`;

	$self->{remove_list_path}="$output_path/length.exclusion.lst";
	
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

