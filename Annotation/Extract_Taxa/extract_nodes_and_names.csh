#!/bin/csh

cut -f 1,3 nodes.dmp > nodes.tsv
cut -f 1,5 nodes.dmp > levels.tsv
grep "scientific name" names.dmp | cut -f 1,3 > names.tsv
