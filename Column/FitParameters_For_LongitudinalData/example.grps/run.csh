#!/bin/csh


../FitParametersForLongitudinalData.r \
	-i Example_PCs.tsv \
        --time_offset_cn=WeekOffset \
        --subject_id_cn=MouseID \
        --target_values_fn=Example.targets \
        --group_cn=SexMet \
        -l -m -d \
        -o Example

