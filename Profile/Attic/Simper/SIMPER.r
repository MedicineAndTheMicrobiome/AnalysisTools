#!/usr/bin/env Rscript

###############################################################################
library(MASS);
library(vegan);
library('getopt');

params=c(
	"summary_table", "i", 1, "character",
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",
	"bootstraps", "b",1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <summary table>\n",
	"	-f <factors>\n",
	"	-b <number of bootstraps>\n",
	"	[-o <output filename root>]\n\n",
	"This  script bootstraps the  samples  present in the summary  file and  calculates  the
 variance  around the percentage  generated  by SIMPER for each  taxon. The  script also
 finds the  Spearman's rank correlation between the taxa abundance and SIMPER percentage.",
	"\n\n");

if(!length(opt$summary_table) || !length(opt$factors) || !length(opt$bootstraps)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputFnameRoot=gsub(".xls", "", opt$summary_table);
}else{
	OutputFnameRoot=opt$outputroot;
}

NumBootstraps=opt$bootstraps;
SummaryTableFname=opt$summary_table;
FactorsFname=opt$factors;


cat("\n");
cat("summary Table Filename: ", SummaryTableFname, "\n", sep="");
cat("Factors Filename: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");

###############################################################################

load_summary_table=function(fname){
	summary_table=as.matrix(read.delim(fname, sep="\t",  header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""));
	mat_dim=dim(summary_table);
	
	return(summary_table);
}

#-------------------------------------------------------------------------------

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

#-------------------------------------------------------------------------------

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE));
	return(factors);
}

##############################################################################

summary_table=load_summary_table(SummaryTableFname);
summary_table=  summary_table[,-1]; #remove the total column
num_samples=nrow(summary_table);
num_taxa=ncol(summary_table);
sample_names=rownames(summary_table);
colnames_summary_table=colnames(summary_table);
sample_totals=apply(summary_table, 1, sum);

cat("Number of rows in summary table=", num_samples,"\n");
cat("Number of columns in summary table=", num_taxa,"\n");

factors=load_factors(FactorsFname);
num_samples_factors=nrow(factors);
#sanity check to make sure number of samples is the same in summary table and in factor file
if(num_samples != num_samples_factors){
	cat("Number of samples in summary table and bootstrap_summary_fctors file do no match\n");
	q(status=-1);
}
factors=factors[,1]; #take only the second column (condition)
factor_levels=levels(factors);
num_factor_levels=length(factor_levels);

normalized_summary_table=normalize(summary_table);

unique=1:num_samples; #for naming
bootstrap_counter= 1;

bootrap_simper_summary_table=matrix(0, nrow=nrow(summary_table), ncol=ncol(summary_table));
bootstrap_simper=matrix(0, nrow=NumBootstraps, ncol=num_taxa); #simper average percentage contribution of each taxa to dissimilarity for every bootstrap

rank_correlation=matrix(0, nrow=1, ncol=NumBootstraps); #rank taxa by simper contribution to dissimilarity
comp_correlation=matrix(0, nrow=1, ncol=NumBootstraps);

###############################################################################

cat("\n\nBootstrapping:\n");
pb <- txtProgressBar(min = 0, max = NumBootstraps, style = 3)

for(bs in 1:NumBootstraps){
	
	setTxtProgressBar(pb, bs);

	unique_counter=1;

	bootstrap_simper_summary_fc=file("temp_bootstrap.xls", "w");
	write(paste("sample_id", paste(colnames(summary_table), collapse="\t"), sep="\t"), file=bootstrap_simper_summary_fc);

	bootstrap_factors_outname=paste(OutputFnameRoot, ".", "bootstrap",NumBootstraps,".factors.txt", sep="");
	bootstrap_factors_fc=file(bootstrap_factors_outname, "w");
	write(paste("sample_id", "condition" , sep="\t"), file=bootstrap_factors_fc)

	# Generate replicates
	samples=sample(1:num_samples, num_samples, replace=TRUE);
	for(i in samples){

		reads_per_sample=sample_totals[i];	

		obs=sample(1:num_taxa, replace=TRUE, reads_per_sample, prob=normalized_summary_table[i,]);
		counts=as.vector(table(c(obs,1:num_taxa))-1);

		name=paste(sample_names[i], ".", unique[unique_counter], ".", bs , sep="");
		unique_counter= unique_counter + 1;

		bootstrap_summary_outline=paste(name, paste(counts, collapse="\t"), sep="\t");
		write(bootstrap_summary_outline, file=bootstrap_simper_summary_fc);

		bootstrap_factors_outline=paste(name, factors[i]);
		write(bootstrap_factors_outline, file=bootstrap_factors_fc);
	}

	close(bootstrap_simper_summary_fc);
	close(bootstrap_factors_fc);

	#perform simper here and calculate get the confidence interval
	bootstrap_simper_summary_table=load_summary_table("temp_bootstrap.xls");
	bootstrap_simper_normalized_summary_table=normalize(bootstrap_simper_summary_table);

	bootstrap_factors=load_factors(bootstrap_factors_outname);
	bootstrap_factors=bootstrap_factors[,1];

	bootstrap_factors_level_assign= list();
	mean_composition_by_factor=list();


	for(i in 1:num_factor_levels){
		bootstrap_factors_level_assign[[i]] = which(bootstrap_factors==factor_levels[i]);#gets index of which group each sample belongs to
		mean_composition_by_factor[[i]] = rep(NaN, num_taxa);#mean for each factor

		for(j in 1:num_taxa){
			mean_composition_by_factor[[i]][j] = mean(bootstrap_simper_normalized_summary_table[bootstrap_factors_level_assign[[i]],j]);#get mean composition for each factor and each taxa
		}
	}

	mean_composition_by_factor_matrix=matrix(unlist(mean_composition_by_factor), ncol=ncol(summary_table), byrow=TRUE);
	mean_composition_by_factor_final=array();
	#get mean composition of entire bootstrap summary table
	for(j in 1:ncol(summary_table)){
			mean_composition_by_factor_final[j] = mean(mean_composition_by_factor_matrix[,j]);#mean for each taxa on bootstrapped data
	}

	ranked_taxa=rank(mean_composition_by_factor_final);#ranked average composition (higher the rank higher the composition)

	sim= simper(bootstrap_simper_normalized_summary_table, bootstrap_factors)

	for(j in 1:(num_factor_levels-1)){
			simper_average=unlist(sim[[1]][2], use.names=FALSE);#average contribution of each taxa to overall dissimilarity
			#simper_order=unlist(sim[[1]][8], use.names=FALSE);#stores the indices of original table in order
			simper_overall=unlist(sim[[1]][3], use.names=FALSE);#overal contribution
			simper_percent=simper_average/simper_overall;
	
			simper_contributions_rank=rank(simper_percent);#ranked average contribution (higher the rank higher the contribution)
	
			bootstrap_simper[bs,]= simper_percent;
	
			#cat("Composition vs Contribution to dissmilarity\n");
			comp_correlation[bs]= cor(mean_composition_by_factor_final, simper_percent, method="pearson");
	
			#cat("Rank abundance vs Rank dissmilarity\n");
			rank_correlation[bs]= cor(ranked_taxa, simper_contributions_rank, method="spearman");
	}
	file.remove(bootstrap_factors_outname);
	#file.remove("temp_bootstrap.xls");
}
close(pb);

sorted_comp_correlation=sort(comp_correlation);
sorted_rank_correlation=sort(rank_correlation);

cat("\nMaking output file...\n\n");

confidence_95=matrix(0, nrow=3, ncol=num_taxa);
sorted_bootstrap_simper=matrix(0, nrow=NumBootstraps, ncol=num_taxa);

#95% confidence 
num_taxa_contributing=1;

for(i in 1:num_taxa){
	sorted_bootstrap_simper[,i]=sort(bootstrap_simper[,i]);

	for(j in 1:3){
		if(j==1){
			confidence_95[j,i]=sorted_bootstrap_simper[ceiling(0.05*NumBootstraps),i];
		}else if(j==2){
			confidence_95[j,i]=sorted_bootstrap_simper[ceiling(0.95*NumBootstraps),i];
		}else if(j==3){
			confidence_95[j,i]=median(sorted_bootstrap_simper[,i]);
		}
	}
}

confidence_interval_outname=paste(OutputFnameRoot, ".", "confidence_interval95.xls", sep="");
confidence_interval_fc=file(confidence_interval_outname, "w");
write(paste("stat", paste(colnames(summary_table), collapse="\t"), sep="\t"), file=confidence_interval_fc);
write(paste("CI95_lower perc dissmilarity", paste(confidence_95[1,], collapse="\t"), sep="\t"), file=confidence_interval_fc);
write(paste("CI95_higher perc dissmilarity", paste(confidence_95[2,], collapse="\t"), sep="\t"), file=confidence_interval_fc);
write(paste("median", paste(confidence_95[3,], collapse="\t"), sep="\t"), file=confidence_interval_fc);
#write(paste("CI95_lower Composition correlation" , comp_correlation[ceiling(0.05*NumBootstraps)] , sep="\t"), file=confidence_interval_fc);
#write(paste("CI95_higher Composition correlation" , comp_correlation[ceiling(0.95*NumBootstraps)] , sep="\t"), file=confidence_interval_fc);
#write(paste("CI95_lower Rank correlation" , sorted_rank_correlation[ceiling(0.025*NumBootstraps)] , sep="\t"), file=confidence_interval_fc);
#write(paste("CI95_higher Rank correlation" , sorted_rank_correlation[ceiling(0.975*NumBootstraps)] , sep="\t"), file=confidence_interval_fc);
close(confidence_interval_fc);

orig_st_mean_composition_by_factor_final=apply(normalized_summary_table, 2, mean);
#print(orig_st_mean_composition_by_factor_final)
orig_st_simper= simper(normalized_summary_table, factors);

orig_st_simper_contributions=array(unlist(orig_st_simper[[1]][2]));#get averages from SIMPER
orig_st_simper_overall=unlist(orig_st_simper[[1]][3]);
#print(orig_st_simper_contributions)
orig_st_simper_percent=orig_st_simper_contributions/orig_st_simper_overall;
orig_st_simper_contributions_rank=rank(orig_st_simper_percent);

orig_st_ranked_taxa=rank(orig_st_mean_composition_by_factor_final);
#print(orig_st_ranked_taxa);

#common_contri_compos=which(orig_st_simper_contributions_rank==orig_st_ranked_taxa);
#print(common_contri_compos);

boxplot_sorted_bootstrap_simper=matrix(0, nrow=NumBootstraps, ncol=20);
boxplot_sorted_bootstrap_simper_col_names=matrix(0,nrow=1, ncol=20)

j=1;
for(i in 1:num_taxa){
		if((orig_st_simper_contributions_rank[i]>(num_taxa-20)) && (j<21)){

		boxplot_sorted_bootstrap_simper[,j]=sorted_bootstrap_simper[,i];
		boxplot_sorted_bootstrap_simper_col_names[j]=colnames_summary_table[i];
		j=j+1;
	}
}

boxplot_sorted_bootstrap_simper_median=apply(boxplot_sorted_bootstrap_simper, 2, median);
#print(boxplot_sorted_bootstrap_simper_median);
boxplot_sorted_bootstrap_simper_order=order(boxplot_sorted_bootstrap_simper_median, decreasing=TRUE);
#print(boxplot_sorted_bootstrap_simper_order);
boxplot_top20=matrix(0, nrow=NumBootstraps, ncol=20);


j=1;
for(i in boxplot_sorted_bootstrap_simper_order){
	boxplot_top20[,j]=boxplot_sorted_bootstrap_simper[,i];
	#barplot_top20_contribution_rank[j]=orig_st_simper_contributions_rank[i];
	#barplot_top20_ranked_taxa[j]=orig_st_ranked_taxa[i];
	j=j+1;
}

orig_st_simper_contributions_rank_order=order(orig_st_simper_contributions_rank, decreasing=TRUE);
orig_st_ranked_taxa_order=order(orig_st_ranked_taxa, decreasing=TRUE);

barplot_top20_contribution_rank=matrix(0, nrow=1, ncol=20);
barplot_top20_contribution_rank_colnames=array(0, dim=20);
barplot_top20_contribution=array(0, dim=20);
barplot_top20_contribution_rank_colnames_reln=list();

j=1;
for(i in orig_st_simper_contributions_rank_order){
	if(j<21){
		barplot_top20_contribution_rank[j]=orig_st_simper_contributions_rank[i];
		barplot_top20_contribution[j]=orig_st_simper_percent[i];
		barplot_top20_contribution_rank_colnames[j]=colnames_summary_table[i];
		barplot_top20_contribution_rank_colnames_reln[[colnames_summary_table[i]]]=orig_st_simper_percent[i];
		j=j+1;
	}
}

barplot_top20_ranked_taxa=matrix(0, nrow=1, ncol=20);
barplot_top20_ranked_taxa_colnames=array(0, dim=20);
barplot_top20_compos_taxa=array(0, dim=20);
#barplot_top20_compos_colnames_reln=list();
j=1;
for(i in orig_st_ranked_taxa_order){
	if(j<21){
		barplot_top20_ranked_taxa[j]=orig_st_ranked_taxa[i];
		barplot_top20_compos_taxa[j]=orig_st_mean_composition_by_factor_final[i];
		barplot_top20_ranked_taxa_colnames[j]=colnames_summary_table[i];
		#barplot_top20_compos_colnames_reln[[colnames_summary_table[i]]]=orig_st_mean_composition_by_factor_final[i];
		j=j+1;
	}
}

barplot_top20_contribution_plot_order=array(0, dim=20);

j=1;
for(i in barplot_top20_ranked_taxa_colnames){
	barplot_top20_contribution_plot_order[j]=barplot_top20_contribution_rank_colnames_reln[[i]];
	j=j+1;
}

barplot_top20_contribution_rank_order=rank(-barplot_top20_contribution_plot_order);

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);
layout(matrix(c(1,2), 2, 1, byrow = TRUE));

boxplot(boxplot_top20, axes=FALSE, ylim=c(0,max(boxplot_top20)+0.05), ylab="Fraction Contribution to Dissimilarity (SIMPER)", main="Top 20 Taxa: Bootstrapped Contribution to Dissimilarity");
axis(2, at=seq(0.00,0.5, by=0.05), labels=seq(0.00,0.5, by=0.05));
#axis(1,at=seq(0, 20, by=1), label=rep("", 11));
axis(1,at=seq(0, 20, by=1), label=c("", boxplot_sorted_bootstrap_simper_col_names), las=2, cex.axis=0.6);

mtext(paste("Spearman's Rank Correlation Coefficient (95% CI) = ","[", round(sorted_rank_correlation[ceiling(0.025*NumBootstraps)]*100, digits=2), ",", round(sorted_rank_correlation[ceiling(0.975*NumBootstraps)]*100, digits=2), "]", sep=" "),side=3,line=0, at=5, cex=0.7);

#legend("topright", boxplot_sorted_bootstrap_simper_col_names, fill=c("springgreen3", "tomato1", "tan2", "steelblue4", "orangered", "lightsalmon1", "chartreuse3", "darkgoldenrod", "darkmagenta", "deeppink4"));

bar_col=rainbow(20);

layout(matrix(c(1,2,3), 3,1, byrow=TRUE));
#par(mfrow=c(3,1));
par(xpd=TRUE);
barplot_taxa_rank=barplot(barplot_top20_compos_taxa, beside=TRUE, ylim=c(0, max(barplot_top20_compos_taxa)+0.07), col=bar_col, ylab="Fraction Composition", main="Top 20 Taxa: Observed Composition");
#text(x=barplot_taxa_rank, y=c(barplot_top20_compos_taxa), labels=c(barplot_top20_ranked_taxa_colnames), cex=0.7, srt=90, offset=0.2, pos=4);
text(x=barplot_taxa_rank, y=c(barplot_top20_compos_taxa), labels=c(num_taxa-barplot_top20_ranked_taxa+1), cex=0.7, pos=3);
par(las=2);
barplot_contribution_rank=barplot(barplot_top20_contribution_plot_order, ylim=c(0, max(barplot_top20_contribution_plot_order)+0.07), ylab="Fraction Contribution to Dissimilarity (SIMPER)",  names.arg=c(barplot_top20_ranked_taxa_colnames), col=bar_col, main="Top 20 Taxa: Observed Contribution to Dissimilarity", cex.names=0.8);
#axis(1,at=seq(0, 20, by=1), label=c("", barplot_top20_ranked_taxa_colnames), las=2, cex.axis=0.6);
#axis(2, at=seq(0.00,0.5, by=0.05), labels=seq(0.00,0.5, by=0.05));
#par(xpd=TRUE);
#text(x=barplot_contribution_rank, y=rep(0,20), labels=c(barplot_top20_ranked_taxa_colnames), cex=0.8, offset=0.2, srt=90, pos=2);
text(x=barplot_contribution_rank, y=c(barplot_top20_contribution_plot_order), labels=c(barplot_top20_contribution_rank_order), cex=0.7, pos=3);

dev.off();