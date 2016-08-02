#!/bin/csh


GeneOntology_toCounts.pl \
	-n ExampleInputFiles/read_to_go_path.example \
	-m ExampleInputFiles/GO.map \
	-c ExampleInputFiles/category_depth \
	-o counts 
