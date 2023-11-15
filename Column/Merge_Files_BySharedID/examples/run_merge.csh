#!/bin/csh

cat <<EOF > targets
16S.min_100.ETA.v1v2.alr.top15.tsv
16S.min_100.GA.v1v2.alr.top15.tsv
16S.min_100.Oral.v1v2.alr.top15.tsv
16S.min_100.RS.v1v2.alr.top15.tsv
EOF

~/git/AnalysisTools/Column/Merge_Files_BySharedID/Merge_Files_BySharedID.r \
	-l targets \
	-o merged.tsv \
	-C SampleID

