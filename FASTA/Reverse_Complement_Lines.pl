#!/usr/bin/env perl

###############################################################################

sub revcomp{

=head2 str revcomp(str sequence);

This function will return the sequence reverse complemented.

=cut
    my $sequence = shift;

    $sequence =~ tr/atugcyrswkmbdhvnxATUGCYRSWKMBDHVNX/taacgryswmkvhdbnxTAACGRYSWMKVHDBNX/;
    my $revcomp = "";
    for (my $j=length($sequence)-1;$j>-1;$j--) {
        $revcomp .= substr($sequence,$j,1);
    }
    return $revcomp;

}

###############################################################################

print STDERR "Enter in nucleotide sequence.  Separate sequences with carriage return.  Control-D to quit.\n";
print STDERR "Note: This does not take a FASTA file.  Single line sequences only.\n";

while(<STDIN>){
	chomp;
	print revcomp($_) . "\n";
}
	
