#!/bin/csh

../Create_Demographics_Table.r \
	-f grads_joined_by_gradid.bal.wCalc.sarc_only.tsv \
	-s smdsc.phengrp
	#-s predniso
	#-s demog.gender
