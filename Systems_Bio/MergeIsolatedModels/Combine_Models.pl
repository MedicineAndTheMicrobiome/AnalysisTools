#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_a $opt_b $opt_o $opt_A $opt_B $opt_e $opt_f $opt_i $opt_x $opt_L);

getopts("a:b:o:A:B:e:f:i:x:L:");

my $UNDEF_LB_LIMIT=0;
my $LB_LIMIT=-1000;
my $UB_LIMIT=1000;

my $usage = "
	usage:
	$0
		-a <Input filename, model A>
		-b <Input filename, model B>	
		-x <Exchange Bounds>

		-o <Output filename root>

		-A <compartment mapping for A>
		-B <compartment mapping for B>

		-e <extracellular for A>
		-f <extracellular for B>
		[-i <ID mapping file>]

		[-L <default lowerbound limit $UNDEF_LB_LIMIT>]

	Generates a combined model.

	1.) Read in reactions
	2.) Map individual exchange reactions into shared space
	3.) Create exchange reactions into shared space

	The ID mapping file will ensure that both models are 
	able to refer to the same metabolites in the shared
	extracellular space.

	Example Invocation:

	Combine_Models.pl \
		-a Rhodoferax_Ferrireducens \
		-b Thiobacillus_Denitrificans \
		-x Exchange_Bounds.tsv \
		-A e=RFe,c=RFc\
		-B e0=TDe,c0=TDc\
		-e RFe \
		-f TDe \
		-i CompoundList.txt.orig \
		-o Combined



";

if(!defined($opt_o) || 
	!defined($opt_a) || !defined($opt_b) || 
	!defined($opt_e) || !defined($opt_f) || 
	!defined($opt_A) || !defined($opt_B)
){
	die $usage;
}

my $orgAfname=$opt_a;
my $orgBfname=$opt_b;
my $orgA_mapping=$opt_A;
my $orgB_mapping=$opt_B;
my $orgA_extracellular=$opt_e;
my $orgB_extracellular=$opt_f;
my $Output=$opt_o;
my $ID_MappingFile=$opt_i;
my $Exchange_Bounds=$opt_x;

if(defined($opt_L)){
	$LB_LIMIT=$opt_L;
}

###############################################################################

print STDERR "Organism A: $orgAfname\n";
print STDERR "Organism B: $orgBfname\n";
print STDERR "Output: $Output\n";
print STDERR "Compartment Mapping A: $orgA_mapping\n";
print STDERR "Compartment Mapping B: $orgB_mapping\n";
print STDERR "Extracellular for A: $orgA_extracellular\n";
print STDERR "Extracellular for B: $orgB_extracellular\n";
print STDERR "ID Mapping File: $ID_MappingFile\n";
print STDERR "\n";

###############################################################################

my $ABBR_COL=0;
my $NAME_COL=1;
my $EQN_COL=2;
my $REV_COL=3;
my $COMP_COL=4;
my $LB_COL=5;
my $UB_COL=6;
my $OBJ_COL=7;
my $RULE_COL=8;
my $SUBS_COL=9;
my $NUM_RXN_COLUMNS=10;
my $SHARED_COMPART_ID="sh";


###############################################################################

my %A_compart_hash;
my %B_compart_hash;

sub parse_compart_mapping{
	my $hash_ref=shift;
	my $mapping_string=shift;

	my @maps=split ",", $mapping_string;
	foreach my $map(@maps){
		my ($src, $dst)=split "=", $map;
		$src=~s/^ +//g;
		$src=~s/ +$//g;
		${$hash_ref}{$src}=$dst;
	}	
}

sub translate_compart{
	my $string_arr_ref=shift;
	my $mapping_ref=shift;
	my $req_brackets=shift;

	my $arr_len=$#{$string_arr_ref}+1;

	for(my $i=0; $i<$arr_len; $i++){
		foreach my $src (keys %{$mapping_ref}){
			my $dst=${$mapping_ref}{$src};
			if($req_brackets){
				${$string_arr_ref}[$i]=~s/\[$src\]/[$dst]/g;
			}else{
				${$string_arr_ref}[$i]=~s/$src/$dst/g;
			}	
		}
	}

}

parse_compart_mapping(\%A_compart_hash, $orgA_mapping);
parse_compart_mapping(\%B_compart_hash, $orgB_mapping);

###############################################################################
# 
# 	GENERATE COMBINED REACTION FILE
#
###############################################################################

sub load_reaction_file{
	my $fname=shift;

	print STDERR "Loading $fname...\n";
	open(REACT_FH, "<$fname") || die "Could not open reaction file for reading: $fname\n";
	
	my @reactions;
	my $rxn_idx;
	while(<REACT_FH>){
		chomp;
		
		my @columns=split "\t", $_;

		# Skip header
		if($columns[$EQN_COL] eq "equation"){
			print STDERR "Header found and removed.\n";
			next;
		}

		# Copy columns into array
		for(my $i=0; $i<=$#columns; $i++){
			$reactions[$rxn_idx][$i]=$columns[$i];
		}

		$rxn_idx++;		
	}

	print STDERR "$rxn_idx reactions loaded.\n";

	return(\@reactions);
	
}

###############################################################################

sub load_id_map{
	my $fname=shift;
	my %idmap;
	
	open(IN_FH, "<$fname") || die "Could not open ID mapping for reading: $fname\n";
	while(<IN_FH>){
		chomp;
		my ($src, $dst)=split "\t";
		if($src eq "" || $dst eq ""){
			next;
		}
		$idmap{$src}=$dst;
	}
	close(IN_FH);
	return(\%idmap);
}

###############################################################################

sub parse_rxn{
	my $in_str=shift;

	my ($global_compart_str, $rxn_str);
	if($in_str=~/ : /){
		($global_compart_str, $rxn_str)=split " : ", $in_str;
	}else{
		$rxn_str=$in_str;
	}
	
	my $dir;
	if($rxn_str=~/-->/){
		$dir="-->";
	}elsif($rxn_str=~/<--/){
		$dir="<--";
	}elsif($rxn_str=~/<==>/){
		$dir="<==>";
	}

	my ($lhs, $rhs)=split $dir, $rxn_str;
	$lhs=~s/^ //;
	$lhs=~s/ $//;
	$rhs=~s/^ //;
	$rhs=~s/ $//;

	my $glob_compart="";
	if($global_compart_str=~/\[(.+)\]/){
		$glob_compart=$1;
	}

	#print STDERR "$lhs, $rhs, $dir, $glob_compart\n";
	return($lhs, $rhs, $dir, $glob_compart);
}

sub is_exchange_reaction{
	my $rhs=shift;
	my $dir=shift;
	if($rhs eq "" && $dir eq "<==>"){
		return(1);
	}else{
		return(0);
	}
}

sub translate_rxn_side{
	my $rxn_str=shift;
	my $shared_id_map_hash=shift;
	my $translated_id_hash=shift;

        foreach my $key(keys %{$shared_id_map_hash}){
                my $value=${$shared_id_map_hash}{$key};
                if(($rxn_str=~s/$key/$value/g)>0){
			print STDERR "$key -> $value\n";
			${$translated_id_hash}{$value}=$key;
		}
        }

        return($rxn_str);
}

sub remap_exchange_and_external_reactions{
	my $rxn_info_ref=shift;
	my $shared_exch_rxn_hash_ref=shift;
	my $shared_id_map_hash=shift;
	my $extracellular_compartment=shift;
	my $num_rxns=$#{$rxn_info_ref}+1;

	my %translated_id_hash;

	for(my $i=0; $i<$num_rxns; $i++){
		my $reaction_str=${$rxn_info_ref}[$i][$EQN_COL];
		my ($lhs, $rhs, $dir, $glob_compart)=parse_rxn($reaction_str);

		if(is_exchange_reaction($rhs, $dir)){

			my $shared_rhs;
			if(defined(${$shared_id_map_hash}{$lhs})){
				$shared_rhs=${$shared_id_map_hash}{$lhs};	
				$translated_id_hash{$shared_rhs}=$lhs;				
			}else{
				$shared_rhs=$lhs;
			}

			# Replace exchange reaction with external-to-shared compartment
			my $new_react="$lhs\[$glob_compart] <==> $shared_rhs\[$SHARED_COMPART_ID]";
			${$rxn_info_ref}[$i][$EQN_COL]=$new_react;
			${$rxn_info_ref}[$i][$COMP_COL].=", $SHARED_COMPART_ID";

			# Control exchanges from shared, so let individual run free
			${$rxn_info_ref}[$i][$LB_COL]=$LB_LIMIT;
			${$rxn_info_ref}[$i][$UB_COL]=$UB_LIMIT;
			#print STDERR "$new_react\n";
			
			# Insert shared compartment exhange reaction
			${$shared_exch_rxn_hash_ref}{"$shared_rhs"}=1	

		}elsif($extracellular_compartment eq $glob_compart){

			$lhs=translate_rxn_side($lhs, $shared_id_map_hash, \%translated_id_hash);			
			$rhs=translate_rxn_side($rhs, $shared_id_map_hash, \%translated_id_hash);
			${$rxn_info_ref}[$i][$EQN_COL]="[$SHARED_COMPART_ID] : $lhs $dir $rhs";
			${$rxn_info_ref}[$i][$COMP_COL]="$SHARED_COMPART_ID";
		}
	}
	
	return(\%translated_id_hash);
}

sub add_shared_reactions{
	my $hash_ref=shift;
	my $rxn_info_ref=shift;
	my $custom_limits_hash_ref=shift;

	my $rxn_info_length=$#{$rxn_info_ref}+1;

	foreach my $met_id(sort keys %{$hash_ref}){

		my ($lb, $ub);
		if(defined(${$custom_limits_hash_ref}{$met_id})){
			($lb, $ub)=split "\t", ${$custom_limits_hash_ref}{$met_id};
		}else{
			($lb, $ub)=($UNDEF_LB_LIMIT, $UB_LIMIT);
			print STDERR "WARNING: Shared flux limits not defined for $met_id (a metabolite from an individual exchange reaction), using $lb to $ub.\n";
		}
	
		my @new_line=(
			"SHARED_EX_$met_id",
			"Shared $met_id exchange",
			"[$SHARED_COMPART_ID] : $met_id <==>",
			"Reversible",
			"$SHARED_COMPART_ID",
			$lb,
			$ub,
			0,
			"",
			""
		);
	
		for(my $i=0; $i<$NUM_RXN_COLUMNS; $i++){
			${$rxn_info_ref}[$rxn_info_length][$i]=$new_line[$i];
		}
		$rxn_info_length++;
		#print STDERR (join "\t", @new_line) . "\n";
	}
}

sub concatenate_reactions{
	my $rxn_info_src_ref=shift;
	my $rxn_info_dst_ref=shift;

	my $src_len=$#{$rxn_info_src_ref}+1;
	my $dst_len=$#{$rxn_info_dst_ref}+1;

	#print STDERR "Src Length: $src_len   Dst Length: $dst_len\n";
	for(my $i=0; $i<$src_len; $i++){
		for(my $j=0; $j<$NUM_RXN_COLUMNS; $j++){
			${$rxn_info_dst_ref}[$dst_len+$i][$j]=${$rxn_info_src_ref}[$i][$j];
		}
	}
}

sub write_reaction_file{
	my $filename=shift;
	my $rxn_info_ref=shift;

	print STDERR "Writing to $filename\n";
	open(FH, ">$filename") || die "Could not open reaction file for writing: $filename\n";

	# Write header
	my @header=("abbreviation", "name", "equation", "reversible", 
		"compartment", "lowbnd", "uppbnd", "obj_coef", "rule", "subsystem");
	print FH (join "\t", @header) . "\n";

	# Write lines
	my $rxn_len=$#{$rxn_info_ref}+1;
	for(my $i=0; $i<$rxn_len; $i++){
		print FH (join "\t", @{${$rxn_info_ref}[$i]}) . "\n";
	}

	close(FH);
	print STDERR "Done writing.\n";

}

sub remap_reaction_compartments{
	my $reaction_info_ref=shift;
	my $compartment_map_hash_ref=shift;

	my $react_info_length=$#{$reaction_info_ref}+1;
	my @reactions;
	my @compartments;

	for(my $i=0; $i<$react_info_length; $i++){
		$reactions[$i]=${$reaction_info_ref}[$i][$EQN_COL];
		$compartments[$i]=${$reaction_info_ref}[$i][$COMP_COL];
	}

	translate_compart(\@reactions, $compartment_map_hash_ref, 1);
	translate_compart(\@compartments, $compartment_map_hash_ref, 0);

	for(my $i=0; $i<$react_info_length; $i++){
		${$reaction_info_ref}[$i][$EQN_COL]=$reactions[$i];
		${$reaction_info_ref}[$i][$COMP_COL]=$compartments[$i];
	}
}

sub rename_reaction_id{
	my $reaction_info_ref=shift;
	my $id=shift;

	my $react_info_length=$#{$reaction_info_ref}+1;

	for(my $i=0; $i<$react_info_length; $i++){
		${$reaction_info_ref}[$i][$ABBR_COL] .= ".$id";
	}
}

sub load_shared_exchange_bounds{
	my $bounds_file=shift;
	my $bounds_hash_ref=shift;

	open(FH, "<$bounds_file") || die "Could not open bounds file for reading: $bounds_file\n";
	while(<FH>){
		chomp;
		my ($met_id, $lb, $ub)=split "\t", $_;
		${$bounds_hash_ref}{$met_id}="$lb\t$ub";
	}
	close(FH);
}

###############################################################################

my $rxn_info_A_ref=load_reaction_file("$orgAfname\_react.tsv");
my $rxn_info_B_ref=load_reaction_file("$orgBfname\_react.tsv");

my $id_map_hash_ref=load_id_map($ID_MappingFile);

my %shared_exchange_bounds_hash;
load_shared_exchange_bounds($Exchange_Bounds, \%shared_exchange_bounds_hash);

###############################################################################

# Build combined reactions file
remap_reaction_compartments($rxn_info_A_ref, \%A_compart_hash);
remap_reaction_compartments($rxn_info_B_ref, \%B_compart_hash);

rename_reaction_id($rxn_info_A_ref, $orgA_extracellular);
rename_reaction_id($rxn_info_B_ref, $orgB_extracellular);


my %shared_exch_rxn_hash;
my $met_A_translated_hash_ref=remap_exchange_and_external_reactions($rxn_info_A_ref, \%shared_exch_rxn_hash, 
	$id_map_hash_ref, $orgA_extracellular);
my $met_B_translated_hash_ref=remap_exchange_and_external_reactions($rxn_info_B_ref, \%shared_exch_rxn_hash, 
	$id_map_hash_ref, $orgB_extracellular);

my @shared_exch_rxns;
add_shared_reactions(\%shared_exch_rxn_hash, \@shared_exch_rxns, \%shared_exchange_bounds_hash);

my @rxn_info_all_ref;
concatenate_reactions($rxn_info_A_ref, \@rxn_info_all_ref);
concatenate_reactions($rxn_info_B_ref, \@rxn_info_all_ref);
concatenate_reactions(\@shared_exch_rxns, \@rxn_info_all_ref);

write_reaction_file("$Output\_react.tsv", \@rxn_info_all_ref);
my $num_reactions=$#rxn_info_all_ref+1;
print STDERR "Num reactions: $num_reactions\n";

###############################################################################
# 
# 	GENERATE COMBINED METABOLITE FILE
#
###############################################################################

sub load_metabolite_file{
	my $fname=shift;
	my $met_hash_ref=shift;

	open(FH, "<$fname") || die "Could not open metabolite file for reading: $fname\n";
	while(<FH>){
		chomp;
		my ($abb, $nam, $comp)=split "\t", $_;
		if($abb eq "abbreviation"){
			next;
		}else{
			if(defined(${$met_hash_ref}{$abb})){
				print STDERR "Metabolite: $abb / \"${$met_hash_ref}{$abb}\" defined more than once (shared) between models. (Probably ok.)\n";
			}else{
				${$met_hash_ref}{$abb}=$nam;
			}
		}
	}
	close(FH);
}

#------------------------------------------------------------------------------
	
sub write_metabolite_file{
	my $fname=shift;
	my $met_hash_ref=shift;

	print STDERR "Writing to $fname\n";
	open(FH, ">$fname") || die "Could not open metabolite file for writing: $fname\n";

	# Write header
	print FH "abbreviation\tname\n";

	# Write lines
	foreach my $id (sort keys %{$met_hash_ref}){
		print FH "$id\t${$met_hash_ref}{$id}\n";
	}
	close(FH);
	print STDERR "Done writing.\n";
}

#------------------------------------------------------------------------------

# Concatentate metabolites file
my %met_hash;
load_metabolite_file("$orgAfname\_met.tsv", \%met_hash);
load_metabolite_file("$orgBfname\_met.tsv", \%met_hash);

# Add met that were translated but whose translated id's were not previously in either models metabolite list
foreach my $trans_hash(($met_A_translated_hash_ref, $met_B_translated_hash_ref)){
	foreach my $key(keys %{$trans_hash}){
		if(!defined($met_hash{$key})){
			$met_hash{$key}=$met_hash{${$trans_hash}{$key}};
		}
	}
}


my $num_mets=keys %met_hash;

write_metabolite_file("$Output\_met.tsv", \%met_hash);
print STDERR "Num metabolites: $num_mets\n";

###############################################################################
# 
# 	GENERATE COMBINED DESCRIPTION FILE
#
###############################################################################

sub load_description_file{
	my $fname=shift;
	my %desc_hash;
	open(FH, "<$fname") || die "Could not open description file for reading: $fname\n";

	$_=<FH>; chomp;
	my @key=split "\t", $_;
	$_=<FH>; chomp;
	my @val=split "\t", $_; 

	for(my $i=0; $i<=$#key; $i++){
		#print STDERR "$key[$i] / $val[$i]\n";
		$desc_hash{$key[$i]}=$val[$i];
	}
	close(FH);
	return(\%desc_hash);	
}

#------------------------------------------------------------------------------

sub write_description_file{
	my $fname=shift;
	my $desc_hash_ref=shift;

	open(FH, ">$fname") || die "Couldn to open $fname for writing.\n";

	print FH (join "\t", (
		"name",
		"id",
		"description",
		"compartment",
		"abbreviation"
	)) . "\n";

	print FH (join "\t", (
		${$desc_hash_ref}{"name"},
		${$desc_hash_ref}{"id"},
		${$desc_hash_ref}{"description"},
		${$desc_hash_ref}{"compartment"},
		${$desc_hash_ref}{"abbreviation"}
	)) . "\n";
	
	close(FH);
}

#------------------------------------------------------------------------------

my $desc_A_hash_ref=load_description_file("$orgAfname\_desc.tsv");
my $desc_B_hash_ref=load_description_file("$orgBfname\_desc.tsv");

my @A_comp=split ", ", ${$desc_A_hash_ref}{"compartment"};
translate_compart(\@A_comp, \%A_compart_hash, 0);
my @B_comp=split ", ", ${$desc_B_hash_ref}{"compartment"};
translate_compart(\@B_comp, \%B_compart_hash, 0);
my @combined_compartments = sort (@A_comp, @B_comp, $SHARED_COMPART_ID);

my @A_abb=split ", ", ${$desc_A_hash_ref}{"abbreviation"};
translate_compart(\@A_abb, \%A_compart_hash, 1);
my @B_abb=split ", ", ${$desc_B_hash_ref}{"abbreviation"};
translate_compart(\@B_abb, \%B_compart_hash, 1);
my @combined_abbreviations= sort (@A_abb, @B_abb, "[$SHARED_COMPART_ID]");



my %combined_hash_ref;
$combined_hash_ref{"name"}=$Output;
$combined_hash_ref{"id"}=$Output;
$combined_hash_ref{"description"}="Combination of $orgAfname and $orgBfname";
$combined_hash_ref{"compartment"}=join ", ", @combined_compartments;
$combined_hash_ref{"abbreviation"}=join ", ", @combined_abbreviations;

write_description_file("$Output\_desc.tsv", \%combined_hash_ref);

###############################################################################
###############################################################################

print STDERR "\n";
print STDERR "Done.\n";








