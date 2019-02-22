#!/bin/csh

if(1)then
../Make_Pair_Mapping_File.r \
	-f np_sinus.crs_obst.first_visit \
	-s SubjectID \
	-t NasalCavity 
endif

if(1)then
../Make_Pair_Mapping_File.r \
	-f first_vs_second_visit.sinus_only \
	-s SubjectID \
	-t Visit 
endif

