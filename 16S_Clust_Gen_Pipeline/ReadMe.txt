
The scripts in this directory can be used in the following order to set up a mothur run.
Some of the steps can be skipped, but for most workflows, they may be necessary.

1.) Merge_Mates.pl

	This script will take a paired fastq file, merge the forward and reverse
	reads to generate a fasta file.

2.) Subsample_FASTAs_in_Directory.pl

	This script will read in a directory full of fasta files, and create another matching
	directory of subsampled fasta files.

	This may be necessary to make a mothur run more tractable.

3.) Generate_SampleID_to_Fasta_Map.pl

	This script will generate a sample id for the set of fasta files in the directory
	based on the fasta file name.

4.) Assign_Sample_IDs_To_Reads.pl

	This script will generate a single monolithic fasta file and a groups file, based on
	the sample ID to fasta file mapping.

5.) Run_Mothur_Steps.pl
	
	This script will automatically run the mothur steps that are necesary and some
	post processing to generate summary tables
