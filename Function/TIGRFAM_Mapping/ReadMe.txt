#!/bin/csh

wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_14.0_INFO.tar.gz
mkdir TIGR_INFO
cd TIGR_INFO
tar -xf ../TIGRFAMs_14.0_INFO.tar.gz
cd ..
cat TIGR_INFO/* > All_TIGR.INFO
