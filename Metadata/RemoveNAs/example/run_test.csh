#!/bin/csh

../Maximize_NonNA.r \
	-i grads_joined_by_gradid.wCalc.byBALKit.tsv \
	-o no_nas.tsv
