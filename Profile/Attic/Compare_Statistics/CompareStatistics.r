#!/usr/local/bin/Rscript

###############################################################################

progname = commandArgs(FALSE)[4]
args = commandArgs(TRUE)

num_args=length(args);
expected_arg_count=2

if(num_args != expected_arg_count){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, 
		"\n\t\t<association file> <data file>",
"

Example Association File:

	Triplet ID,Lesion,Normal,Control
	1,PL040-01,PN040-01,CE047-02
	2,PL052-01,PN052-01,CB053-02
	3,PL032-01,PN032-01,CA019-01
	4,PL027-01,PN027-01,CD023-01
	7,PL053-01,PN053-01,CB017-01

Example Attribute Data File:

	Sample ID,Shannon,Simpson,Evenness
	PN051-01,2.42149268656223,0.86260510944723,0.71195300208624
	PN007-01,2.40955111407414,0.879917927531631,0.739557924682419
	PN027-01,2.40303752591784,0.871188818169162,0.789298674939225
	PL005-01,2.38907469062623,0.87881060791298,0.811385363232274

This program will compute the change in attributes between the different samples grouped in the 
association file.  
"
	);

	writeLines(usage)

	if(num_args>expected_arg_count){
		writeLines("Too many arguments. ");
	}else{
		writeLines("Not enough arguments. ");
	}

	writeLines(paste("Only", expected_arg_count, "arguments expected.\n"));

	quit(status=0)
}

###############################################################################
# Main program loop

association_file=args[1];
attribute_file=args[2];
output_file=paste(attribute_file, ".pdf", sep="");

cat("\n");
cat("Association File: ", association_file, "\n");
cat("Attribute File:   ", attribute_file, "\n");
cat("Output File:      ", output_file, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load association data

association_matrix<-as.matrix(read.table(association_file, sep=",", header=TRUE, check.names=FALSE))
#print(association_matrix);
num_associations=nrow(association_matrix);
num_association_types=ncol(association_matrix)-1;
cat("Num associations: ", num_associations, "\n");
cat("Num association types: ", num_association_types, "\n");

association_list=vector("pairlist",0);
for(association_id in 1:num_associations){
	association_list[[association_matrix[association_id,1]]]=association_matrix[association_id,2:(num_association_types+1)];
}
#print(association_list);
association_names=names(association_matrix[1,2:ncol(association_matrix)]);

###############################################################################

# Load sample attribute data

attributes_matrix<-as.matrix(read.table(attribute_file, sep=",", header=TRUE, check.names=FALSE))
#print(attributes_matrix);

# Convert attributes into a "hash"/list
num_samples=nrow(attributes_matrix);
num_attributes=ncol(attributes_matrix)-1;
attribute_names=names(attributes_matrix[1,1:num_attributes+1]);
cat("Num samples: ", num_samples, "\n");
cat("Num attributes: ", num_attributes, "\n");
#print(attribute_names);

sample_attribute_list=vector("pairlist",0);
for(sample_id in 1:num_samples){
	sample_attribute_list[[attributes_matrix[sample_id,1]]]=as.real(attributes_matrix[sample_id,2:(num_attributes+1)]);
}
#print(sample_attribute_list);

###############################################################################
# Compute the association type comparisons we need to do, eg. A vs B, B vs C, A vs C, etc.

comp_idx=1;
assoc_comparisons=vector("pairlist",0);
for(i in 1:num_association_types){
	for(j in i:num_association_types){
		if(i!=j){
			assoc_comparisons[[comp_idx]]=c(i,j);
			comp_idx=comp_idx+1;
		}
	}
}
num_comparisons=length(assoc_comparisons);
#print(assoc_comparisons);

###############################################################################
# Cycle through the associations
assoc_ids=names(association_list);

pdf(output_file,width=11,height=8.5)
par(mfrow=c(3,3));

for(attribute in 1:num_attributes){
	cat("Working on attribute: ", attribute_names[attribute] , "\n", sep="");

	all_of_attribute=vector("pairlist",0);
	all_of_compare_str=vector("pairlist",0);
	all_of_normality=vector("pairlist",0);

	for(comparison_idx in 1:num_comparisons){
		comparison=assoc_comparisons[[comparison_idx]];
		compare_str=paste(association_names[comparison[1]], " - ", association_names[comparison[2]], 
			" [", attribute_names[attribute], "]\n");
		#cat("Working on: ", compare_str, "\n", sep="");

		delta=rep(0,1);
		num_deltas=1;
		attrib1=rep(0,1);
		attrib2=rep(0,1);
		for(assoc_id in 1:num_associations){
			#cat("    Working on: ", assoc_ids[assoc_id], " [", assoc_id, "] \n");

			# pick the patient group/pairing
			cur_assocs=association_list[[assoc_id]];
			sampleid1=(cur_assocs[comparison[1]]);
			sampleid2=(cur_assocs[comparison[2]]);
			#cat(sampleid1, " vs ", sampleid2, "\n");

			# Get the attributes associate with each patient
			sample_attributes1=sample_attribute_list[[sampleid1]];
			sample_attributes2=sample_attribute_list[[sampleid2]];

			# Get the attribute we want to compare
			test_attribute1=sample_attributes1[attribute];
			test_attribute2=sample_attributes2[attribute];
			
			# Do the math and store those with valid/non-Null attributes
			if(length(test_attribute1)==1 && length(test_attribute2)==1){
				attrib1[num_deltas]=test_attribute1;
				attrib2[num_deltas]=test_attribute2;
				delta[num_deltas]=test_attribute1-test_attribute2;
				num_deltas=num_deltas+1;
			}else{
				cat(assoc_ids[assoc_id], " Rejected.\n");
			}	
		}

		# Store computed differences so we can plot then together later.
		all_of_attribute[[comparison_idx]]=delta;
		all_of_compare_str[[comparison_idx]]=compare_str;
		shap1=shapiro.test(attrib1);
		shap2=shapiro.test(attrib2);
		all_of_normality[[comparison_idx]]=c(shap1$p, shap2$p);

		# Plot overlayed density functions
		dens_mix=density(c(attrib1,attrib2));
		dens1=density(attrib1, bw=dens_mix$bw);
		dens2=density(attrib2, bw=dens_mix$bw);
		all_hist=hist(c(attrib1,attrib2), plot=FALSE);
		max=max(c(dens1$y,dens2$y));
		plot(dens1, col="blue", type="l", main=attribute_names[attribute], ylim=c(0,max), xlab="");
		lines(dens2, col="red" );
		mtext( paste(association_names[comparison[1]]), col="blue", side=3, line=-1, cex=.7);
		mtext( paste(association_names[comparison[2]]), col="red", side=3, line=-2, cex=.7);

		# Plot side-by-side histograms
		bothhist=hist(c(attrib1,attrib2),plot=FALSE)
		ahist=hist(attrib1,breaks=bothhist$breaks,plot=FALSE)
		bhist=hist(attrib2,breaks=bothhist$breaks,plot=FALSE)
		bothmat=t(matrix(c(ahist$intensities, bhist$intensities),length(bothhist$breaks)-1,2))
		dimnames(bothmat)=list(NULL, bothhist$breaks[1:(length(bothhist$breaks)-1)])
		barplot(bothmat,beside=TRUE,col=c("blue","red"), ylim=c(0,max(bothmat)*1.1),space=c(0,.3), main=attribute_names[attribute])
		mtext( paste(association_names[comparison[1]]), col="blue", side=3, line=-1, cex=.7);
		mtext( paste(association_names[comparison[2]]), col="red", side=3, line=-2, cex=.7);

	}

	all_differences=numeric(0);
	all_max=0;
	for(comparison_idx in 1:num_comparisons){
		all_differences=c(all_differences,all_of_attribute[[comparison_idx]]);
		tmp_hist=hist(all_of_attribute[[comparison_idx]], plot=FALSE, breaks=10);
		max_count=max(tmp_hist$density);
		if(all_max<max_count){
			all_max=max_count;
		}
	}
	x_range=c(min(all_differences), max(all_differences));
	x_diff=abs(x_range[1]-x_range[2]);
	x_range=c(x_range[1]-.2*x_diff, x_range[2]+.2*x_diff);


	for(comparison_idx in 1:num_comparisons){
		delta=all_of_attribute[[comparison_idx]];
		title=all_of_compare_str[[comparison_idx]]

		hist(delta, main=title, ylim=c(0,all_max*1.2), xlim=x_range, freq=FALSE,breaks=10);

		# For pair t-test
		shap=shapiro.test(delta);
		d=mean(delta);
		s=sd(delta);
		n=length(delta);
		t=d/(s/sqrt(n));
		p=1-pt(abs(t),n-1);

		# For wilcoxon signed rank test
		w=wilcox.test(delta);

		lines(density(delta));
		mtext(
			paste(
				"Mean = ", sprintf("%4.4f", d),
				"  ",
				"Stdev = ", sprintf("%4.4f", s),
				"  ",
				"N = ", sprintf("%i", n)
			), side=3, line=1, cex=.7);

		norm=all_of_normality[[comparison_idx]];
		mtext(paste("Shapiro-Wilks p-value (normality): ", sprintf("%4.2f",norm[1]), "/", sprintf("%4.2f",norm[2])), side=3, line=0, cex=.7);
		mtext(
			paste(
				"Paired t-test: t = ", sprintf("%4.4f", t),
				"  ",
				"p-value = ", sprintf("%4.4f", p)
			), side=3, line=-1, cex=.7);

		mtext(
			paste(
				"Wilcoxon Signed Rank Test: p-value = ", sprintf("%4.4f", w$p.value)
			), side=3, line=-2, cex=.7);


		cat("\n");
	}
}

dev.off();

warnings_str=warnings();
if(!is.null(warnings_str)){
	print(warnings());
}

###############################################################################

writeLines("Done.\n")

q(status=0)
