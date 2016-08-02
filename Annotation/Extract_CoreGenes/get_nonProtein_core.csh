#!/bin/csh

cut -f 1,4,6 Core_TIGRFAM_GO_EC.tsv | grep -v "Protein synthesis" | grep -v "Protein fate" | cut -f 1,3 > _nonProtein.tsv
cut -f 1 _nonProtein.tsv > nonProt.hmms.tsv
grep -f nonProt.hmms.tsv Core_Roles.tsv > _nonProteinIds

cut -f 1 _nonProteinIds | sort -u > nonProt.genes.tsv
cut -f 3 _nonProteinIds | sort -u > nonProt.go.tsv
cut -f 2 _nonProtein.tsv | grep -v NULL | sort -u > nonProt.ecs.tsv

rm _nonProteinIds
rm _nonProtein.tsv
