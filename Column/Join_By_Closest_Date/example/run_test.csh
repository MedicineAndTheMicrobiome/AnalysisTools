#!/bin/csh

~/git/AnalysisTools/Column/Join_By_Closest_Date/Join_By_Closest_Date.r \
	--main_filen=main_data.tsv \
        --main_keycn=subject_id \
        --main_offcn=collection_date_form \
        --main_dateformat="%Y-%m-%d" \
        --aux_filen=aux_data.tsv \
        --aux_keycn=subject_id \
        --aux_offcn=Date \
        --aux_dateformat="%Y-%m-%d" \
        --aux_valuecn=BMI \
        --outputfn=combined.tsv

