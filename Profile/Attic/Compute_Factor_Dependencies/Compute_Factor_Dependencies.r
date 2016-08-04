#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
options(useFancyQuotes=F);

params=c(
	"factors", "f", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors>\n",
	"\n",
	"This script will read in the factors and try to calculate either\n",
	"dependency between factors, or correlation between factors.\n",
	"\n",
	"If both factors are categorical, then there will be a test of dependence.\n",
	"If both factors are ordinal, then there will be a test of correlation.\n",
	"If a factor is ordinal and the other is categorical, then the ordinal\n",
	"	factor will be binned, and then a test of dependence will be made.\n",
	"\n");

if(!length(opt$factors)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factors;

OutputRoot=FactorsFile;
OutputRoot=gsub("\\.txt","", OutputRoot, ignore.case = T)
OutputRoot=gsub("\\.tsv","", OutputRoot, ignore.case = T)

cat("\n");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE));
	return(factors);
}

plot_correl_heatmap=function(mat, title=""){

	par(family="Courier");
	par(oma=c(10, 10, 1, .5));
	par(mar=c(5.1, 4.1, .5, .5));

	# Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

	# Remember that rows and columsn are reversed in the image
	image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

	# Pad strings
	cnames=paste(colnames(mat), " ", sep="");
	rnames=paste(rownames(mat), " ", sep="");
	
	# Get longest length of each column or row label
	cname_max_len=max(nchar(cnames));
	rname_max_len=max(nchar(rnames));

	# Get the number of rows and columns
	ncols=ncol(mat);
	nrows=nrow(mat);

	cscale=min(c(25/cname_max_len, 25/ncols));
	rscale=min(c(25/rname_max_len, 25/nrows));

        max_width=max(nchar(sprintf("%.2f",mat)));
	cell_cex=(3.5/max_width)*sqrt(min(c(cscale, rscale))^2);

        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){
			str=sprintf("%.2f",mat[i,j]);
			str=gsub("0\\.",".", str);
			used_cell_cex=cell_cex;
			if(nchar(str)>3){
				used_cell_cex=cell_cex*4/nchar(str);
			}
                        text(i,j,labels=str, cex=used_cell_cex, srt=45);
                }
        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

	# Plot the title
	mtext(title, line=0, at=nrows*.5, side=3, font=2);

}


plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
	par(mar=rep(0,4));

	num_lines=length(strings);
	
	top=max(as.integer(num_lines), 52);

	plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);
	for(i in 1:num_lines){
		#cat(strings[i], "\n", sep="");
		text(0, top-i, strings[i], pos=4, cex=.8); 
	}
}

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
num_factors=ncol(factors);
num_samples=nrow(factors);

cat("Num Factors: ", num_factors, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");

summary(factors);

factor_names=colnames(factors);
factor_sample_names=rownames(factors);

	
##############################################################################

pdf(paste(OutputRoot, ".dependencies.pdf", sep=""), height=11, width=8.5);

# Plot Histograms for all factors, even if they are categories
par(mfrow=c(4,3));

categorical=rep(T, num_factors);
for(i in 1:num_factors){
	
	#cat("Factor Name: ", factor_names[i], "\n");
	val=(factors[,i]);

	# Determine if factor is categorical or numerical
	att=attributes(val);
	if(is.null(att)){
		categorical[i]=F;
	}

	t=table(val);

	if(categorical[i]){
		t_sort=sort(t, decreasing=T);

		# Compute shrinkage for many samples
		labsize_h=min(1, 25/length(t_sort));

		# Compute shrinkage for  long label names
		max_label_length=max(nchar(names(t_sort)));
		labsize_l=min(1, 10/max_label_length);		
		lmar=(5.1*max_label_length/6)
		lmar=min(lmar, 15);

		labsize=min(c(labsize_h, labsize_l));
		
		# Increase bottom margin based on largest label size
		orig_mar=par()$mar;
		temp_mar=orig_mar;
		temp_mar[1]=lmar;
		par(mar=temp_mar);

		# Generate barplot
		barplot(t_sort, main=factor_names[i], xlab="", ylab="Frequency", 
			col="white", las=2, cex.names=labsize);

		# Return to default margins
		par(mar=orig_mar);
	}else{
		h=hist(val, main=factor_names[i], xlab="");

		num_cat=length(t);
		if(num_cat < 5){
			cat_names=names(t);
			tot=sum(t);
			norm=t/tot;
			for(i in 1:num_cat){
				mtext(
					paste(cat_names[i], sprintf(": %3.2f%%", 100*norm[i]), sep=""),
					line=-i,
					cex=.75
				);
			}
		}
	}
}

##############################################################################

pvalues=matrix(NA, ncol=num_factors, nrow=num_factors);
correl=matrix(NA, ncol=num_factors, nrow=num_factors);

colnames(pvalues)=factor_names;
rownames(pvalues)=factor_names;

colnames(correl)=factor_names;
rownames(correl)=factor_names;

num_tests=(num_factors*(num_factors-1))/2;
cat("Num Tests Performed: ", num_tests, "\n");
pvalue_vect=numeric(num_tests);
k=0;

for(i in 1:num_factors){
	for(j in i:num_factors){

		# Only compute half matrix, copy the other half

		cat(factor_names[i], " vs ", factor_names[j], "\n", sep=""); 

		fA=categorical[i];
		fB=categorical[j];

		valA=factors[,i];
		valB=factors[,j];

		if(!fA && !fB){

			cat("\tCorrelation:\n");
			rho=cor(valA, valB);
			res=cor.test(valA, valB);

			pvalues[i,j]=res$p.value;
			correl[i,j]=rho;

			#cat(rho, " / ", res$p.value, "\n");

		}else if(fA && fB){

			cat("\tDependence:\n");
			t=table(valA, valB);
			res=chisq.test(t, simulate.p.value=T);
			pvalues[i,j]=res$p.value;

			#print(t);
			#print(res);
			
		}else{
			cat("\tDiscretized Dependence:\n");

			if(!fB){
				temp=valB;
				valB=valA;
				valA=temp;			
			}

			h=hist(valA, plot=F);

			whichmid=numeric(num_samples);	
			for(ix in 1:num_samples){
				whichmid[ix]=max(which(h$breaks<=valA[ix]));
				if(is.na(whichmid[ix])){
					whichmid[ix]=length(h$mids);
				}
				#cat("Value: ", valA[i], " Bin: ", h$mid[whichmid[i]], "\n");	
			}

			t=table(whichmid, valB);
			res=chisq.test(t, simulate.p.value=T);
			pvalues[i,j]=res$p.value;

			print(t);
			print(res);
		}

		pvalues[j,i]=pvalues[i,j];
		correl[j,i]=correl[i,j];

		k=k+1;
		pvalue_vect[k]=pvalues[i,j];

	}

}

##############################################################################

par(mfrow=c(1,1));
plot_correl_heatmap(correl, title="Pearson Correlation Coefficients");
plot_correl_heatmap(pvalues, title="P-values: Null Hypothesis is No correlation/dependence");
plot_correl_heatmap(log(pvalues, 10), title="Log10(p-values)");

# Adjust p-vlaues for multiple test
adjusted_pvalue_vect=p.adjust(pvalue_vect, "holm");
#print(cbind(pvalue_vect, adjusted_pvalue_vect));

# Place corrected pvalues back into 2D matrix
adjusted_pvalue_matrix=matrix(NA,ncol=num_factors, nrow=num_factors);
colnames(adjusted_pvalue_matrix)=factor_names;
rownames(adjusted_pvalue_matrix)=factor_names;
k=1;
for(i in 1:num_factors){
	for(j in i:num_factors){
		adjusted_pvalue_matrix[i,j]=adjusted_pvalue_vect[k];		
		adjusted_pvalue_matrix[j,i]=adjusted_pvalue_matrix[i,j]
		k=k+1;
	}
}

# Plot heatmap for BH corrected pvalues
plot_correl_heatmap(adjusted_pvalue_matrix, title="Holm Corrected P-values");

# Plot heatmap for significant correlations
significant_matrix=adjusted_pvalue_matrix>.05;
diag(significant_matrix)=1;
plot_correl_heatmap(significant_matrix, title="Significant Correl/Depend after Holm Correction (p<.05)");

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
