#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_e $opt_d);
getopts("ed");

###############################################################################

print STDERR "

If you're waiting, and nothing is happening, try piping a file in.

This program will automatically convert quality values from 0 to 60, to 
ASCII characters from '?' to '{', or vice versa.  It will automatically 
detect which conversion to do.  

The following options are available:
	-e Expect quality values to encode
	-d Expect ASCII values to decode

These are useful to check for errors.  Otherwise, errors will be silent.

";

###############################################################################

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
my $i=0;
while(<STDIN>){
	chomp;
	
	if(/^>/){
		print STDOUT "$_\n";
	}else{
		process_record($_);
	}

	if(!($i%50000)){
		print STDERR ".";
	}
	$i++;
}
close(FASTA_FH);

print STDERR "\nCompleted.\n";

###############################################################################

sub process_record{
	my $line=shift;
	my $outstr="";
	my $asci;
	my $qual;

	if($line=~/\s+/){
		#encode
		if(defined($opt_d)){
			die "You expected decoding, but these look like they are already decoded.\n";
		}

		my @qualvals=split /\s+/, $line;
		foreach $qual(@qualvals){
			if($qual ne ""){

				if($qual>60){
					die "Your quality value has exceeded the logical max of 60.\n";
				}elsif($qual<0){
					die "Your quality value has exceeded the logical min of 0.\n";
				}

				$asci=$qual+ord('?');
	
				$outstr.=chr($asci);
			}
		}
	}else{
		#decode
		if(defined($opt_e)){
			die "You expected encoding, but these look like they are already encoded.\n";
		}

		my $strlen=length($line);
		for(my $i=0; $i<$strlen; $i++){
			$asci=ord(substr($line, $i, 1));

			$qual=$asci-ord('?');

			if($qual>60){
				die "Your quality value has exceeded the logical max of 60.\n";
			}elsif($qual<0){
				die "Your quality value has exceeded the logical min of 0.\n";
			}

			$outstr.= sprintf("%02i ", $qual);
		}
	}
	print STDOUT "$outstr\n";
}

###############################################################################

