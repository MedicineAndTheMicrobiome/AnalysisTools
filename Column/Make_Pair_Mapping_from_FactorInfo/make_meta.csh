#!/bin/csh

~/git/AnalysisTools/Column/Make_Pair_Mapping_from_FactorInfo/Paired_to_Metadata.r \
	-p example/np_sinus.crs_obst.first_visit.NasalCavity.map.tsv \
	-o example/np_sinus.crs_obst.first_visit.NasalCavity.meta.tsv 

~/git/AnalysisTools/Column/Make_Pair_Mapping_from_FactorInfo/Paired_to_Metadata.r \
	-p example/np_sinus.crs_obst.first_visit.NasalCavity.map.wNAs.tsv \
	-o example/np_sinus.crs_obst.first_visit.NasalCavity.wNAs.meta.tsv

~/git/AnalysisTools/Column/Make_Pair_Mapping_from_FactorInfo/Paired_to_Metadata.r \
	-p example/np_sinus.crs_obst.first_visit.NasalCavity.map.tsv \
	-o example/np_sinus.crs_obst.first_visit.NasalCavity.meta.location.tsv \
	-c location
