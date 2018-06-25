#!/bin/csh

R CMD build Bookmarker

R CMD REMOVE Bookmarker
R CMD INSTALL Bookmarker_0.1.tar.gz
