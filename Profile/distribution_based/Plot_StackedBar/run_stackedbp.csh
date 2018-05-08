#!/bin/csh

set summary_table=#../../PACTLHIVCKOW_OW.taxa.genus.no_mcl.summary_table.tsv
set factor_file=#../../pact_lhiv_carl.wform
set out=$summary_table:t:r:r

~/git/AnalysisTools/Profile/distribution_based/Plot_StackedBar/Plot_StackedBar.r \
        -i $summary_table \
        -f $factor_file \
        -o $out \
        -s ";" \

       #-t <top categories to display, default=4>
       #-l <label abundances greater than specified threshold, default=1.0, recommended 0.01>
       #-c <crossing/interactions list, e.g., "var1,var2,var3" >
       #-d <diversity, default=tail.>

