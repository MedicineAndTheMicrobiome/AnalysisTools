#!/usr/bin/env perl

###############################################################################

while(<STDIN>){

	chomp;
	my ($otu, $count, $taxa)=split "\t", $_;

	if($otu eq "OTU"){
		print "$_\n";
	}else{
	
		#print "$taxa\n";
		my @breakdown=split ";", $taxa;

		my @new_breakdown;
		foreach my $component(@breakdown){
			#print "\t$component\n";	
			
			$component=~s/"//g;
			my $id;
			my $perc;
			if($component ne "unclassified"){
				if($component=~/(\S+)\((\d+)\)/){
					$id=$1;
					$perc=$2;
				}else{
					print STDERR "Error parsing: $otu / $count / $taxa\n";
				}
			
				if($id ne "unclassified"){
					push @new_breakdown, $id;
				}
			}
		}

		my $new_taxa_name=join " ", @new_breakdown;

		#print "'$new_taxa_name'\n";
		#print "\n\n";

		my $new_out=join "\t", ($otu, $count, $new_taxa_name);

		print "$new_out\n";
	}

}

