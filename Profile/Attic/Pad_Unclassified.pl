#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_b $opt_o);
use FindBin;
use FileHandle;
getopts("f:b:o:");

my $TAXA_LIST_EXT="taxa_list";
my $SUMMARY_TABLE_XLS_EXT="summary_table.xls";
my $R_INPUT_TXT_EXT="r_input.txt";
my $DEFAULT_BERGEY="$FindBin::Bin/bergeyTrainingTree.xml";

my $usage = "usage:
$0
	[-f <_countcum.txt file>]
	[-b <bergeyTrainingTree.xml>]
	[-o <padded countcum.txt file>]

	if -f not specified, input read from STDIN.
	if -o not specified, output sent to STDOUT.

	Default BergeyFile: $DEFAULT_BERGEY
";

# If a parent is inserted, then its name is taken from the child
# Example parent insertion:
#
#	Original:
#		3       Bacteria OD1
#		3       Bacteria OD1 OD1_genera_incertae_sedis
#
#	Inserted:
#		3       Bacteria OD1 family:OD1_genera_incertae_sedis#uncertain_suborder
#		3       Bacteria OD1 genus:OD1_genera_incertae_sedis#uncertain_subfamily
#		3       Bacteria OD1 order:OD1_genera_incertae_sedis#uncertain_subclass
#		3       Bacteria OD1 subclass:OD1_genera_incertae_sedis#uncertain_class
#		3       Bacteria OD1 subfamily:OD1_genera_incertae_sedis#uncertain_family
#		3       Bacteria OD1 suborder:OD1_genera_incertae_sedis#uncertain_order
#
# If a child is inserted, then its name is inherited from it's parent.
# Example children insertion:
#
#	Original:
#		2000	Bacteria
#		894	Bacteria Fusobacteria
#		894     Bacteria Fusobacteria Fusobacteria
#		894     Bacteria Fusobacteria Fusobacteria Fusobacteriales
#		894     Bacteria Fusobacteria Fusobacteria Fusobacteriales Fusobacteriaceae
#		894     Bacteria Fusobacteria Fusobacteria Fusobacteriales Fusobacteriaceae Fusobacterium
#
#	Inserted:
#		1106    Bacteria domain:Bacteria#uncertain_class
#		1106    Bacteria domain:Bacteria#uncertain_family
#		1106    Bacteria domain:Bacteria#uncertain_genus
#		1106    Bacteria domain:Bacteria#uncertain_order
#		1106    Bacteria domain:Bacteria#uncertain_phylum
#		1106    Bacteria domain:Bacteria#uncertain_subclass
#		1106    Bacteria domain:Bacteria#uncertain_subfamily
#		1106    Bacteria domain:Bacteria#uncertain_suborder
#
# The colon (:), should be read as \"of\".

print STDERR $usage;

my $input_count_cumu_file=$opt_f;
my $output_count_cumu_file=$opt_o;
my $bergey_class_file;

if(!defined($opt_b)){
	$bergey_class_file=$DEFAULT_BERGEY;
}else{
	$bergey_class_file=$opt_b;
}

# Load tree structure from bergey tree definition
my ($id_to_rank_ref, $parent_to_children_ref, $child_to_parent_ref, $id_to_name_ref, $rank_to_id_ref)=load_bergey_xml_file($bergey_class_file);

my %id_to_rank=%{$id_to_rank_ref};
my %parent_to_children=%{$parent_to_children_ref};
my %child_to_parent=%{$child_to_parent_ref};
my %id_to_partialname=%{$id_to_name_ref};
my %rank_to_id=%{$rank_to_id_ref};

# Generate id-to-fullname and fullname-to-it lookups
my ($fullname_to_id_ref, $id_to_fullname_ref)=generate_fullname_to_taxaid_map($child_to_parent_ref, $id_to_name_ref, $id_to_rank_ref);
my %fullname_to_id=%{$fullname_to_id_ref};
my %id_to_fullname=%{$id_to_fullname_ref};

# Load counts (referenced by full name)
my $count_hash_ref=load_count_cumu($input_count_cumu_file, $fullname_to_id_ref);
my %count_hash=%{$count_hash_ref};

#------------------------------------------------------------------------------

# Define the order of the ranks
my @ranks=("domain", "phylum", "class", "subclass", "order", "suborder", "family", "subfamily", "genus");
my %rank_order=(
	"domain"=>0,
	"phylum"=>1,
	"class"=>2,
	"subclass"=>3,
	"order"=>4,
	"suborder"=>5,
	"family"=>6,
	"subfamily"=>7,
	"genus"=>8
);

#----------------------------------------------------------------------------------

print STDERR "\n\nFixing counts that need to be interpolated bottom up.  \n  (eg.  genus assigned without an upper rank assigned)\n";
# Fix bottom up missing counts
# In this bottom up code, we need to make modifications to the tree, but for the top down, we can make inferences about the counts
#  without modifying the tree.
for(my $rank_idx=$#ranks; $rank_idx>0; $rank_idx--){

	print STDERR "\t$ranks[$rank_idx]\n";
	foreach my $child_id(@{$rank_to_id{$ranks[$rank_idx]}}){
		
		# Get parent of node
		my $parent_id=$child_to_parent{$child_id};

		# Compute the rank that the parent *should* be.
		my $immediate_parent_rank_numerical=$rank_order{$id_to_rank{$child_id}}-1;
		my $immediate_parent_rank_text=$ranks[$immediate_parent_rank_numerical];

		# If the parent skips a rank, then we need to insert a new parent and rebuild links.
		if($rank_order{$id_to_rank{$parent_id}} != $immediate_parent_rank_numerical){

			# Make new unique parent id
			my $new_parent_id=$child_id . ":$immediate_parent_rank_text";

			# Grab the partial name with the finest granularity
			my $child_fullname=$id_to_fullname{$child_id};
			my @child_name_comp=split / /, $child_fullname;
			my $partial_name=pop @child_name_comp;

			if($partial_name=~/(.+):(.+)#(.+)/){
				$partial_name=$2;
			}

			# Make new parent partial name:  <parent rank>:<child partial name>#uncertain_<parent rank>
			my $new_parent_partial_name=$id_to_rank{$child_id} . ":" . $partial_name . "#uncertain_$immediate_parent_rank_text";
			push @child_name_comp, $new_parent_partial_name;
			my $new_parent_fullname=join " ", @child_name_comp;
			$id_to_fullname{$new_parent_id}=$new_parent_fullname;

			# INSERT NEW PARENT NODE INTO TREE
			# Point new parent at child's former parent
			$child_to_parent{$new_parent_id}=$child_to_parent{$child_id};
			# Point new parent at child
			push @{$parent_to_children{$new_parent_id}}, $child_id;
			# Point child's former parent to new parent
			push @{$parent_to_children{$parent_id}}, $new_parent_id;
			# Point child at new parent
			$child_to_parent{$child_id}=$new_parent_id;

			# Assign rank to new parent
			$id_to_rank{$new_parent_id}=$immediate_parent_rank_text;
			# Save parent's partial name
			$id_to_partialname{$new_parent_id}=$new_parent_partial_name;
			# Save parent's full name
			$id_to_fullname{$new_parent_id}=$new_parent_fullname;
			# Save parent's id in list of ids for it's rank
			push @{$rank_to_id{$immediate_parent_rank_text}}, $new_parent_id;

			# Bump up counts from child to parent
			$count_hash{$new_parent_id}=$count_hash{$child_id};
		}
	}
}
print STDERR "\n";

#----------------------------------------------------------------------------------

print STDERR "Removing children that are not exactly one rank down.\n\n";
# Delete children that are not one rank away. When we inserted new parents, we didn't
# have the former parents stop pointing at the former child.
for(my $rank_idx=0; $rank_idx<$#ranks; $rank_idx++){

	my $target_rank=$ranks[$rank_idx];
	foreach my $parent_id(@{$rank_to_id{$target_rank}}){

		my $parent_rank_numerical = $rank_idx;

		# Keep track of which chilren we want to keep
		my @immediate_children;
		foreach my $child_id(@{$parent_to_children{$parent_id}}){
			my $child_rank_numerical = $rank_order{$id_to_rank{$child_id}};
			if($child_rank_numerical - $parent_rank_numerical == 1){
				push @immediate_children, $child_id;
			}#else{
			# Child is abandoned.
			#}
		}

		# Reset list of children for this parent
		$parent_to_children{$parent_id}=\@immediate_children;
	}
}

#----------------------------------------------------------------------------------

# If a rank has an undefined count, take the sum of it's children

for(my $rank_idx=$#ranks; $rank_idx>0; $rank_idx--){

	my $target_rank=$ranks[$rank_idx];

	print STDERR "Rank: $target_rank\n";

	foreach my $parent_id(@{$rank_to_id{$target_rank}}){

		#print STDERR "   Looking at $parent_id / $id_to_fullname{$parent_id}\n";
		my $parent_count=$count_hash{$parent_id};
		if(!defined($parent_count) || $parent_count==0){



			my $sum=0;
			foreach my $child_id(@{$parent_to_children{$parent_id}}){
				$sum+=$count_hash{$child_id};
			}

			if($sum>0){
				print STDERR "    Updating sum: $sum\n";
				$count_hash{$parent_id}=$sum;
			}
		}
	}
}


#----------------------------------------------------------------------------------

print STDERR "\n\nFixing counts that can be interpolated top down.\n  (eg. where finer granularity classification was not possible.)\n\n";

# Fix top down missing counts
for(my $rank_idx=0; $rank_idx<$#ranks; $rank_idx++){

	my $target_rank_text=$ranks[$rank_idx];
	my $below_rank=$ranks[$rank_idx+1];

	print STDERR "\t$target_rank_text\n";
	
	# For every parent, if the sum of the children's count is less than the parent, then insert a new child to hold the difference.
	foreach my $parent_id (@{$rank_to_id{$target_rank_text}}){

		# Get parent's count
		my $parent_count=$count_hash{$parent_id};

		# If the parent has a count then check it's children.
		if($parent_count>0){

			# Sum up the counts of the children
			my $child_sum=0;
			my @children=@{$parent_to_children{$parent_id}};
			foreach my $child_id (@children){
				# Only sum up children that are directly/immediately below the parent in rank
				if($id_to_rank{$child_id} eq $below_rank){
					my $child_count=$count_hash{$child_id};
					$child_sum+=$child_count;
				}
			}

			# If the sum of the children is less than parent's, adopt a placeholder child
			if($child_sum<$parent_count){

				# Compute counts that should be in the adopted child
				my $diff=$parent_count-$child_sum;

				# Assign name to adopted child based on parent's fullname
				my $prior_fullname=$id_to_fullname{$parent_id};

				# Actually, we aren't inserting a child, we are just inserting new entries into the count hash
				# For every placeholder children, insert decendents all the way down.
				for(my $insIdx=$rank_idx+1; $insIdx<=$#ranks; $insIdx++){
					# Assign name: <parent full name> <actual rank>:<partial name>#uncertain_<last known rank>
 					my $new_fullname="$prior_fullname $target_rank_text:$id_to_partialname{$parent_id}#uncertain_$ranks[$insIdx]";
					my $new_id="$parent_id:$ranks[$insIdx]";
					# Store new count
					$id_to_fullname{$new_id}=$new_fullname;
					$count_hash{$new_id}=$diff;
				}
			}
		}
	}

}

#--------------------------------------------------------------------------------
# Output the results

# Pick a destination to write the counts
my $fh=new FileHandle;
if(!defined($output_count_cumu_file)){
	$fh=\*STDOUT;
}else{
	$fh->open(">$output_count_cumu_file") || die "Could not open $output_count_cumu_file for reading.\n";	
}

# Go through every name and print out the counts
foreach my $id (sort keys %{count_hash}){
	# Only report non zero counts
	if($count_hash{$id}>0){
		print {$fh} "$count_hash{$id}\t$id_to_fullname{$id}\n";
	}
}

print STDERR "Done.\n\n";

###############################################################################

sub generate_fullname_to_taxaid_map{

# Generate the full name by climbing up the bergey tree for every node

	my $child_to_parent_hash_ref=shift;
	my $id_to_name_hash_ref=shift;
	my $id_to_rank_hash_ref=shift;

	my %child_to_parent_hash=%{$child_to_parent_hash_ref};
	my %id_to_name_hash=%{$id_to_name_hash_ref};
	my %id_to_rank_hash=%{$id_to_rank_hash_ref};

	my %fullname_to_taxaid_hash;
	my %subfreefullname_to_taxaid_hash;
	my %taxaid_to_fullname_hash;

	foreach my $child(keys %child_to_parent_hash){

		my @names;
		my @subfree_names;
				
		my $taxid=$child;
		do{
			my $name=$id_to_name_hash{$taxid};
			unshift @names, $name;
			
			if(!($id_to_rank{$taxid}=~/^sub/)){
				unshift @subfree_names, $name;
			}


			$taxid=$child_to_parent_hash{$taxid};
		}while($taxid!=0);

		my $whole_name_string=join " ", @names;
		my $whole_subfree_name_string=join " ", @subfree_names;

		$fullname_to_taxaid_hash{$whole_name_string}=$child;
		$taxaid_to_fullname_hash{$child}=$whole_name_string;

		# Include a mapping from subclass free names to id.
		# This is because the input may be missing subclasses

		if(!($id_to_rank_hash{$child}=~/^sub/)){
			$fullname_to_taxaid_hash{$whole_subfree_name_string}=$child;
		}
	}

	return(\%fullname_to_taxaid_hash, \%taxaid_to_fullname_hash, \%subfreefullname_to_taxaid_hash);

}

###############################################################################

sub load_count_cumu{

# Load the count_cum file into a fullname -> count hash

	my $path=shift;
	my $fullname_to_id_ref=shift;
	my %fullname_to_id=%{$fullname_to_id_ref};

	my %count_hash;

	my $fh=new FileHandle;
	if(!defined($path)){
		$fh=\*STDIN;
	}else{
		$fh->open("<$path") || die "Could not open $path for reading.\n";	
	}

	while(<$fh>){
		chomp;

		$_=~s/"//g;
		my @line=split /\t/, $_;

		my $total=shift @line;
		my $leaf=pop @line;
		if($leaf ne "(any)"){
			push @line, $leaf;
		}

		my $taxa=join " ", @line;

		if(defined($fullname_to_id{$taxa})){
			my $id=$fullname_to_id{$taxa};
			$count_hash{$id}=$total;
		}else{
			die "Could not identify taxa: '$taxa'\n";
		}
	}

	return(\%count_hash);
}

###############################################################################

sub load_bergey_xml_file{

# Load the tree structure from the bergey xml

	my $bergey_xml_file=shift;
	my %id_to_rank;
	my %parent_to_children;
	my %child_to_parent;
	my %id_to_name;
	my %rank_to_id;

	open(FH, "<$bergey_xml_file") || die "Could not open $bergey_xml_file";
	while(<FH>){

		chomp;
		$_=~s/^\<TreeNode //;
		$_=~s/\>\<\/TreeNode\>$//;

		if($_=~/name="(.+)" taxid="(\d+)" rank="(\S+)" parentTaxid="(\d+)" leaveCount="(\d+)" genusIndex="([\d-]+)"/){
			my ($name, $taxid, $rank, $parent_taxid, $leave_count, $genus_index)=($1, $2, $3, $4, $5, $6);

			# Clean up invalid characters
			$name=~s/&quot;//g;
			$name=~s/ /_/g;

			# Store all taxa for a particular rank 
			$id_to_rank{$taxid}=$rank;

			# Store parent to child relationship
			push @{$parent_to_children{$parent_taxid}}, $taxid;

			# Store parent for every child		
			$child_to_parent{$taxid}=$parent_taxid;

			# Store name for every taxid
			$id_to_name{$taxid}=$name;

			# Store rank to taxid
			push @{$rank_to_id{$rank}}, $taxid;

		}else{
			if(!(($_=~/trainsetNo/) || ($_=~/name="Root"/))){
				print STDERR "$_\n";
				print STDERR "\tparse error\n";
			}
		}
	}
	close(FH);

	# Insert Unknown into tree
	$id_to_rank{"Unknown"}="domain";
	$id_to_name{"Unknown"}="Unknown";
	$child_to_parent{"Unknown"}=0;
	push @{$rank_to_id{"domain"}}, "Unknown";

	return(\%id_to_rank, \%parent_to_children, \%child_to_parent, \%id_to_name, \%rank_to_id);
}


###############################################################################



