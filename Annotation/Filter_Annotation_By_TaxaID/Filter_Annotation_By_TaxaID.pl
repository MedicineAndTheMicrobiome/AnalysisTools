#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_t $opt_a $opt_f $opt_o);

getopts("t:a:f:o:");

my $usage = "
	usage:
	$0

	-t <taxa class file>
	-a <annotation file>
	-f <taxa id's to filter>
	-o <output filtered annotation filename root>

	Taxa class file:
		query_id, taxa_id,
                superkingdom, phylum, class, order, family, genus, species

	Annotation file:
		query_id, ...<remaining annotation>


	Taxa id's to filter:
		taxa_id1\\ttaxa_name1\\n
		taxa_id1\\ttaxa_name2\\n
		...
		taxa_idn\\ttaxa_namen\\n

	There will be two filtered annotation files generated.  One containing
	the reads that were filtered, and the other containing those kept.


";

###############################################################################

if(!defined($opt_t) || 
	!defined($opt_a) ||
	!defined($opt_f) ||
	!defined($opt_o)
){
	die $usage;
}

my $TaxaClassFile=$opt_t;
my $AnnotationFile=$opt_a;
my $TaxaFilterFile=$opt_f;
my $OutputFilenameRoot=$opt_o;

###############################################################################

sub load_hash{

	my $fname=shift;
	my %hash;
	
	open(FH, "<$fname") || die "Could not open $fname\n";

	while(<FH>){
		chomp;
		my ($item, $comment)=split "\t", $_;
		$hash{$item}=1;
	}
	
	close(FH);

	return(\%hash);	
}

###############################################################################

my $SUPERKINGDOM_COL=2;
my $SPECIES_COL=9;

sub filter_load_taxa_class{
	
	my $taxa_classes_fname=shift;
	my $taxa_hash_ref=shift;
	my %filter_hash;

	open(FH, "<$taxa_classes_fname") || die "Could not open $taxa_classes_fname\n";

        while(<FH>){

                chomp;
                my @taxa_ids=split "\t", $_;

		my $filt_query=0;

		for(my $i=$SUPERKINGDOM_COL; $i<=$SPECIES_COL; $i++){
			
			if(${$taxa_hash_ref}{$taxa_ids[$i]}==1){
				$filt_query=1;
				last;
			}
		}

		if($filt_query){
			#print "$taxa_ids[0]\n";
			$filter_hash{$taxa_ids[0]}=1;
		}

        }

        close(FH);
	return(\%filter_hash);
}

###############################################################################

# Load taxa id's to filter
my $taxa_id_filt_hash_ref=load_hash($TaxaFilterFile);
print STDERR "Taxa ID's to Filter:\n";
foreach my $taxa_id(keys %{$taxa_id_filt_hash_ref}){
	print STDERR "'$taxa_id'\n";
}
print STDERR "\n";

# Generate hash of reads to filter
my $remove_hash_ref=filter_load_taxa_class($TaxaClassFile, $taxa_id_filt_hash_ref);


# Go through annoation file and split out which annotations to keep
my $outkeepfname="$OutputFilenameRoot\.kept";
my $outfiltfname="$OutputFilenameRoot\.remo";

open(FH, "<$AnnotationFile") || die "Could not open $AnnotationFile\n";
open(KEEP_FH, ">$outkeepfname") || die "Could not open $outkeepfname.\n";
open(FILT_FH, ">$outfiltfname") || die "Could not open $outfiltfname.\n";

my $kept=0;
my $filt=0;
while(<FH>){
	my $orig=$_;
	chomp;
	my @fields=split "\t", $_;
	
	#print("$fields[0]\n");
	if(defined(${$remove_hash_ref}{$fields[0]})){
		print FILT_FH $orig;
		$filt++;
	}else{
		print KEEP_FH $orig;
		$kept++;
	}

}

close(KEEP_FH);
close(FILT_FH);

###############################################################################

print STDERR "Kept: $kept\n";
print STDERR "Removed: $filt\n";
print STDERR "Total: ", ($kept+$filt), "\n";

###############################################################################

print STDERR "Done.\n";

