#!/bin/csh

# These should be in your .cshrc 
setenv PATH ${PATH}:/usr/local/common
source /usr/local/common/env/sybase.csh
setenv SYBASE /usr/local/packages/sybase
setenv DSQUERY SYBPROD

# This is the actual command to run
cat roles.core.sql | runsql -D common -P access -l access > Core_Roles.tsv
cat TIGRFAM_GO_EC.core.sql | runsql -D common -P access -l access > Core_TIGRFAM_GO_EC.tsv

# Replace - with more specific ECs:
#
mv Core_TIGRFAM_GO_EC.tsv Core_TIGRFAM_GO_EC.orig.tsv
cat Core_TIGRFAM_GO_EC.orig.tsv | \
    sed "s/3.1.25.-/3.1.25.1/" | \
    sed "s/2.7.7.-/2.7.7.6/" | \
    sed "s/6.1.1.-/6.1.1.12/" | \
    sed "s/6.3.4.-/6.3.4.19/" | \
    sed "s/3.6.5.-/3.6.5.n1/" > \
	Core_TIGRFAM_GO_EC.tsv

