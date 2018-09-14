#!/bin/csh

dos2unix example.mod

./split_model.pl \
	-m example.mod \
	-d .

