
COMPARING THE TAXONOMIC DIVERSITY BETWEEN MULTIPLE SAMPLES

The purpose of the code in this directory is to visualize the similarities 
or differences between a collection of samples.  This production code was 
adapted from Shibu's original set of shell scripts so that it could be 
easily run on any kind of comparison of samples.

The code runs relatively quickly, since there is no sequence comparison 
involved.  The analysis is performed on the number of cumulative members 
at a particular taxonomic rank, between samples.  These numbers are a 
result of the site specific diversity analysis (Aaron Halpern's Pipeline).  
The input is the *_countcum.txt files.

The first script to run is:

	Pull_Taxonomic_Crosssection_By_Rank.pl

This will produce a list of taxas at a particular rank, a summary table of 
counts, and a simplified table for R's consumption.  The only complicated 
part of this script is building the "sample id"-to-"file name" mapping 
file.  The format is described in the usage of the script.

The next script to run is:

	Graph_Taxonomic_Diversity.r

This will read in the simplified table and produce a dendrogram and a MDS 
plot.

You can find additional information when you run either of the scripts on 
the command line without arguments.

Note that the order of the how the samples are listed in the first script 
will affect the MDS and dendrogram plots.  There's not much you can do 
about it since it is the nature of the algorithm, but be aware if you are 
running any regression tests that order will affect the results subtly.
