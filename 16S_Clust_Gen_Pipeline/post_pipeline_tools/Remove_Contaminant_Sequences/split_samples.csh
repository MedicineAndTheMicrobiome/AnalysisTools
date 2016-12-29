#!/bin/csh

grep -vE "ExtrNegCtrl|GsulfrdCP2|Banthracis|PCRNegCtrl|BCTL" sample_ids > expm_samples.list
grep -E "ExtrNegCtrl|GsulfrdCP2|Banthracis|PCRNegCtrl|BCTL" sample_ids > ctrl_samples.list
