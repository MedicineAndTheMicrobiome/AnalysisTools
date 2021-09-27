#!/bin/csh

../FitParametersForLongitudinalData.r \
	-i BMI_over_time.tsv \
	--time_offset_cn=Days \
	--subject_id_cn=SubjectID \
	--target_values_fn=target_list \
	-l \
	-m \
	-d \
	-r \
	-n \
	-o BMI_over_time.param.tsv 
