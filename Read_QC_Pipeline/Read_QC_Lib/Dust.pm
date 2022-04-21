#!/usr/bin/env perl

###############################################################

package Read_QC_Lib::Dust;

=pod

=head1 NAME

Read_QC_Lib::Dust

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
	$self->{output_directory}="";
	$self->{processed_fasta}="";
	$self->{executable_path}="";
	$self->{stderr_log}="stderr";
	$self->{stdout_log}="stdout";
	$self->{diagnostic_exec_path}="";
	$self->{remove_list_path}="";

	# Dust specific
	$self->{input_format}="fasta";		# Could be fastq
	$self->{output_format}="remove list";	# Could be fastq
	$self->{cutoff}="";
	$self->{DUST_OUT_NAME}="dust.out.fasta";
	$self->{STATS_OUT_NAME}="dust.out.stats";
	$self->{REMOVE_LIST_NAME}="dust.remove.list";
	$self->{purge}=0;
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

sub get_remove_list_path{
	my $self = shift;
	return($self->{remove_list_path});
}

sub set_cutoff{
	my $self = shift;
	$self->{cutoff}=shift;
}

sub get_output_format{
        my $self = shift;
        return($self->{output_format});
}

sub get_input_format{
        my $self = shift;
        return($self->{input_format});
}

sub set_purge{
	my $self = shift;
	$self->{purge}=shift;
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

	# Confirm necessary variables are set
	my $err=0;
	$err+=check_variable($exec_path, "Execution Path");
	$err+=check_variable($output_path, "Result Path");
	$err+=check_variable($fasta_path, "Input FASTA Path");
	die "Variables undefined." unless !$err;

	# Make output directory if necessary
	if(!(-e $output_path)){
		mkdir $output_path;
	}

	# Construct execution string
	my $execute_analysis_string="$exec_path -in $fasta_path -out $output_path/$self->{DUST_OUT_NAME} -outfmt fasta";

	# Execute
	print STDERR "Executing: $execute_analysis_string\n";
	my $res=`$execute_analysis_string`;
	
	return;
}	

sub parse_results{
	my $self=shift;

	my $err=0;
	$err+=check_variable($self->{output_path}, "Result Path");
	die "Variables undefined." unless !$err;

	my $output_file="$self->{output_path}/$self->{DUST_OUT_NAME}";
	my $stats_file="$self->{output_path}/$self->{STATS_OUT_NAME}";

	print STDERR "Reading: $output_file\n";
	print STDERR "Writing: $stats_file\n";
	open(FH, "<$output_file") || die "Could not open $output_file for parsing.\n";
	open(STATS_FH, ">$stats_file") || die "Could not open $stats_file for writing.\n";

	# Define what to do to each fasta record
	my $cnt=0;
	sub process_record{
		my $defline = shift;
		my $sequence = shift;

		$defline=~/^>(\S+)/;
		my $id=$1;
		
		my $seq_len=length($sequence);
		my $num_Ns=($sequence=~tr/a-z//);

		my $perc_hq_seq;
		if($seq_len>0){
			$perc_hq_seq=($seq_len-$num_Ns)/$seq_len;
		}else{
			$perc_hq_seq=0;
		}
		
		my $outline=join "\t", ($id, $num_Ns, $seq_len, sprintf("%3.3f", $perc_hq_seq));
		print STATS_FH "$outline\n";
		$cnt++;
	}

	# Parse file as FASTA
	my ($defline, $prev_defline, $sequence);
	while(<FH>){
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

	close(STATS_FH);
	close(FH);
	print STDERR "Records parsed: $cnt\n";

	return;
}

sub identify_remove_list{
	my $self=shift;

	my $err=0;
	$err+=check_variable($self->{output_path}, "Result Path");
	$err+=check_variable($self->{cutoff}, "Dust Cutoff");
	die "Variables undefined." unless !$err;

	die "Dust Cutoff must be betwen 0 and 1, not $self->{cutoff}" unless $self->{cutoff}<1;

	my $stats_file="$self->{output_path}/$self->{STATS_OUT_NAME}";
	my $cutoff=$self->{cutoff};
	my $remove_list_path="$self->{output_path}/$self->{REMOVE_LIST_NAME}";


	open(STATS_FH, "<$stats_file") || die "Could not open $stats_file\n";
	open(RM_LIST_FH, ">$remove_list_path") || die "Could not open $remove_list_path\n";

	my $num_rem=0;
	while(<STATS_FH>){
		chomp;
		my ($defline, $num_Ns, $seq_len, $perc_hq_seq)=split "\t", $_;
		if($perc_hq_seq<$cutoff){
			print RM_LIST_FH "$defline\n";
			$num_rem++;
		}	
	}

	close(RM_LIST_FH);
	close(STATS_FH);

	print STDERR "Num suggested for removal: $num_rem\n";
	$self->{remove_list_path}=$remove_list_path;
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

	# 2. Parse results into stats
	print STDERR "\nParsing Results:\n";
	$self->parse_results();

	# 3. Identify losers
	print STDERR "\nIdentifying Remove List:\n";
	$self->identify_remove_list();

	# 4. Optionally, generate stats
	print STDERR "\nRunning Diagnostics:\n";
	$self->execute_diagnostics();

	# 5. Purge intermediate files
	if($self->{purge}){
		print STDERR "Purging Dust temporary fasta file: $self->{output_path}/$self->{DUST_OUT_NAME}\n";
		unlink "$self->{output_path}/$self->{DUST_OUT_NAME}" || die "Could not purge $self->{DUST_OUT_NAME}\n";
	}else{
		print STDERR "Leaving temporary file: $self->{output_path}/$self->{DUST_OUT_NAME}\n";
	}	


	print STDERR "\nDone.\n";

	return;
}

###############################################################################

sub DESTROY {
    my $self = shift;
}

###############################################################################

1;

