#!/bin/csh

Extract_GO_Info.pl -m example_input/gene_ontology_ext.obo -o GO.map

cut -f 1,3 GO.map > GO_Info.tsv
