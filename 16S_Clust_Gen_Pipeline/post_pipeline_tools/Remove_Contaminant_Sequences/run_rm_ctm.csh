#!/bin/csh

set PROJDIR = Results_shortnames_750

./Remove_Contaminant_Sequences.pl \
	-a ${PROJDIR}/grads.unique.good.filter.unique.precluster.pick.an.list \
	-g ${PROJDIR}/grads.groups \
	-n ${PROJDIR}/grads.unique.good.filter.unique.precluster.pick.names \
	-c ./ctrl_samples.list \
	-f ${PROJDIR}/grads.fasta \
	-o ./grads_ctmf 
