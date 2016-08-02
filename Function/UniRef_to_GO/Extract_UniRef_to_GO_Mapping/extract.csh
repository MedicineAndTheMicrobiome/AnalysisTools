#!/bin/csh

UniRef_to_Columns.pl -m example_input/idmapping_selected.tab.TEST -o uniref_to_go.map.TEST

UniRef_to_Columns.pl -m example_input/idmapping_selected.tab -o uniref_to_go.map.full

