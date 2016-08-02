#!/bin/csh

set db = ExampleDatabaseFiles

# Extract out: Reaction ID, EC, Enzyme IDs, Participating Pathways IDs
MetacycDAT_to_Columns.pl -m $db/reactions.dat -f "UNIQUE-ID,EC-NUMBER,ENZYMATIC-REACTION,IN-PATHWAY" -o reaction_info.tsv

# Extract out: Enzyme IDs to Common Name
MetacycDAT_to_Columns.pl -m $db/enzrxns.dat -f "UNIQUE-ID,COMMON-NAME" -o enzyme_commonname.tsv

# Extract out: Pathway ID, Pathway Type, Pathway Common Name
MetacycDAT_to_Columns.pl -m $db/pathways.dat -f "UNIQUE-ID,TYPES,COMMON-NAME" -o pathway_type_commonname.tsv

# Extract out: Classes, Class Types, Class Common Name
# Most of these are unused.  But this a tree can be constructed from these relationships.
MetacycDAT_to_Columns.pl -m $db/classes.dat -f "UNIQUE-ID,TYPES,COMMON-NAME" -o object_classes.tsv
