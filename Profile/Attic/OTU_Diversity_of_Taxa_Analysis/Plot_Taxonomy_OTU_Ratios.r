#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2014 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library(plotrix);
library('getopt');

params=c(
	"otu_map", "m", 1, "character",
	"output_filename_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-m <OTU to Taxa Mapping File>\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will go through the top taxa and generate\n",
	"a plot that indicates the distribution of taxa and generate\n",
	"a table that summarizes the diversity of OTUs with each taxa.\n",
	"\n",
	"An example input file looks like:\n",
"
OTU     Size    Taxonomy
Otu0001 1470    Bacteria Proteobacteria Gammaproteobacteria Alteromonadales Alteromonadaceae Alishewanella
Otu0002 625     Bacteria
Otu0003 554     Bacteria Proteobacteria Betaproteobacteria Rhodocyclales Rhodocyclaceae
Otu0004 457     Bacteria Gemmatimonadetes Gemmatimonadetes Gemmatimonadales Gemmatimonadaceae Gemmatimonas
Otu0005 401     Bacteria
Otu0006 389     Bacteria Actinobacteria Actinobacteria Solirubrobacterales
Otu0007 387     Bacteria Cyanobacteria_Chloroplast Chloroplast Chloroplast
Otu0008 351     Bacteria Acidobacteria Acidobacteria_Gp1 Granulicella
Otu0009 301     Bacteria Cyanobacteria_Chloroplast Chloroplast Chloroplast
Otu0010 286     Bacteria Actinobacteria Actinobacteria Actinomycetales Pseudonocardiaceae Actinomycetospora
Otu0011 230     Bacteria
Otu0012 196     Bacteria Bacteroidetes Sphingobacteria Sphingobacteriales Cytophagaceae Spirosoma
Otu0013 174     Bacteria Bacteroidetes Bacteroidia Bacteroidales Porphyromonadaceae
Otu0014 136     Bacteria
Otu0015 119     Bacteria Actinobacteria Actinobacteria
Otu0016 117     Bacteria
Otu0017 107     Bacteria Bacteroidetes Sphingobacteria Sphingobacteriales Cytophagaceae
Otu0018 106     Bacteria
Otu0019 77      Bacteria Firmicutes Clostridia Clostridiales Lachnospiraceae
...
\n",
	"\n",
	"\n");

if(!length(opt$otu_map)){
	cat(usage);
	q(status=-1);
}

MapFile=opt$otu_map;
OutputFilenameRoot=opt$output_filename_root;

cat("Input Map File: ", MapFile, "\n", sep="");
cat("Output Filename Root: ", OutputFilenameRoot, "\n", sep="");

###############################################################################

load_map_file=function(fname){
	mat=as.matrix(read.table(fname, header=T, sep="\t"));
	return(mat);
}

###############################################################################

plot_counts=function(counts_mat, num_taxa=30, num_otus=30){

	# counts_mat[otus, taxa]
	
	par(mar=c(1, 30, 1, 1));
	plot(0,0, type="n", xlim=c(0,num_otus), ylim=c(-num_taxa, 0), yaxt="n", ylab="");

	# Draw circles of area proportional to the OTUs percent contribution to the taxa
	scale=30/num_taxa
	for(i in 1:num_taxa){
		tot=sum(counts_mat[,i]);
		for(j in 1:num_otus){
			radius=sqrt((counts_mat[j, i]/tot)/pi);
			if(radius>0){
				draw.circle(j, -i, radius*scale);
			}

		}
	}

	taxa_names=colnames(counts_mat);
	axis(side=2, at=-1*(1:num_taxa), labels=taxa_names[1:num_taxa], las=2, cex.axis=.7);

}

###############################################################################

compute_taxa_diversity=function(counts_mat){
	
	num_taxa=ncol(counts_mat);
	max_otus=nrow(counts_mat);

	entropy=function(p){
		p=p[p>0];
		ent=-sum(p*log(p));
		return(ent);
	}

	# Create empty matrix
	info_matrix=matrix(0, nrow=num_taxa, ncol=6);
	rownames(info_matrix)=colnames(counts_mat);
	colnames(info_matrix)=c("Entropy", "ProbOfMax", "HasDominant", "NumOTUs", "TotalReads", "MedianReadsPerOTU");

	# Analyze each taxa individually
	for(i in 1:num_taxa){
		#print(counts_mat[,i]);
		x=counts_mat[,i];
		prob=x/(sum(x));
		ent=entropy(prob);
		dom=prob[1]>=.5;
		num_otus=sum(prob>0);
		prob_max=prob[1];
		tot_reads=sum(x);
		med_rpo=median(x[x>0]);
		info_matrix[i,]=c(ent, prob_max, dom, num_otus, tot_reads, med_rpo);
	}

	#print(info_matrix);

	return(info_matrix);

}

###############################################################################

mat=load_map_file(MapFile);

counts=as.numeric(mat[,2]);
taxa_names=mat[,3];

#print(counts);
#print(taxa_names);

unique_taxa_names=unique(taxa_names);
num_uniq_taxa_names=length(unique_taxa_names);

#print(unique_taxa_names);

# Group counts by taxa name;
counts_by_taxa=list();
totals_by_taxa=rep(0, num_uniq_taxa_names);
names(totals_by_taxa)=unique_taxa_names;
max_otus_per_taxa=0;

# Group OTUs by taxa
for(taxa in unique_taxa_names){
	taxa_counts=counts[taxa_names == taxa];
	num_otus_per_taxa=length(taxa_counts);
	max_otus_per_taxa=max(max_otus_per_taxa, num_otus_per_taxa);
	counts_by_taxa[[taxa]]=sort(taxa_counts, decreasing=T);
	totals_by_taxa[taxa]=sum(taxa_counts);
}

cat("Max different OTUs per taxa: ", max_otus_per_taxa, "\n");

# Sort taxa by counts
totals_by_taxa=sort(totals_by_taxa, decreasin=T);
sorted_names=names(totals_by_taxa);
#print(totals_by_taxa);

# Convert into square matrix sorted by most abundant taxa
counts_matrix=matrix(0, nrow=max_otus_per_taxa, ncol=num_uniq_taxa_names);
colnames(counts_matrix)=sorted_names;
for(i in 1:num_uniq_taxa_names){
	cur_name=sorted_names[i];
	cur_counts=counts_by_taxa[[cur_name]];
	counts_matrix[1:length(cur_counts), cur_name]=cur_counts;	
}

# Generate Plots
pdf(paste(OutputFilenameRoot, ".taxa_otu.pdf", sep=""), height=11, width=8.5);
plot_counts(counts_matrix, num_taxa=30, num_otus=15);
plot_counts(counts_matrix, num_taxa=50, num_otus=10);
plot_counts(counts_matrix, num_taxa=100, num_otus=5);

# Compute taxa that have the most diversity and have a dominant OTU
diversity_table=compute_taxa_diversity(counts_matrix);

# Output statistics on each taxa
write.table(diversity_table, file=paste(paste(OutputFilenameRoot, ".taxa_otu.tsv", sep="")),
	sep="\t", col.names=NA, row.names=T);

###############################################################################

dev.off();
cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}
q(status=0)
