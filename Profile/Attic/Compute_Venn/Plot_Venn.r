#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=2

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table1.xls> <Input summary_table2.xls> [<Input summary_table3.xls>]\n\n",
		"This script will sum up all the counts across samples in each summary table and then \n",
		"apply percentage cutoffs across all taxonomic classification in order to generate multiple venn diagrams.",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined (You need to specify at least 2 files).\n")
	quit(status=0)
}

###############################################################################
# Main program loop

cutoffs=c(0,.01,.02,.03,.04,.05);
cutoffs_str=sprintf("%0.2f", cutoffs);
print(cutoffs_str);

file_list=list();
filenames=character(0);
loaded_data=list();
all_categories_list=list();

library(plotrix);

plot_venn = function(count_list, names, cutoff){

	#print(count_list);
	#print(names);
	num_sets=length(names);

	offset=0;
	if(num_sets==3){
		offset=-.5;
	}

	plot(0,0, 
		xlim=c(-2.5,2.5), ylim=c(-2.5-offset,2.5-offset), 
		col="white",  ylab="", xlab="", xaxt="n", yaxt="n",
		bty="n"
		);

	# Draw circles
	# A
	draw.circle(-.6,0,1, col=rgb(1,0,0,.3));
	# B
	draw.circle(.6,0,1, col=rgb(0,0,1,.3));

	# Label counts 
	text(-1,offset*1.8, count_list$A);
	text(1,offset*1.8, count_list$B);
	text(-2,offset*1.8, names[1], srt=-45);
	text(2,offset*1.8, names[2], srt=45);
	text(0,offset*1.8, count_list$AB);

	if(num_sets==3){
		draw.circle(0,1.35,1, col=rgb(0,1,0,.3));
		text(0,3.1, names[3]);
		text(0,2.4, count_list$C);
		text(0,.3, count_list$ABC);
		text(.6,.8, count_list$BC);
		text(-.6,.8, count_list$AC);
	}

	# Cutoff
	mtext(paste("Cutoff = ", cutoff), at=-2.5, side=3, font=2);
}
	
arg_count=1;
while(!(is.na(args[arg_count]))){

	InputFileName=args[arg_count];
	cat("\n   Input File Name: ", InputFileName, "\n")

	###############################################################################
	###############################################################################

	# Load data
	inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))
	#cat("Original Matrix:\n")
	#print(inmat);
	InputFileName=basename(InputFileName)
	# Get input dimensions
	num_orig_samples=nrow(inmat);
	num_orig_categories=ncol(inmat)-1;
	cat("Num samples:", num_orig_samples, "\n");
	cat("Num categories:", num_orig_categories, "\n");

	# Extract out counts, ignore totals
	count_mat=inmat[,2:ncol(inmat)];

	# Get sample/category names
	orig_sample_names=rownames(count_mat);
	#print(orig_sample_names);
	orig_category_names=colnames(count_mat);
	#print(orig_category_names);

	# Only keep the taxon name with the finest classification
	shortened_category_names=rep("",num_orig_categories);
	for(catidx in 1:num_orig_categories){
		#print(category_names[catidx]);
		components=unlist(strsplit(orig_category_names[catidx], " "));
		shortened_category_names[catidx]=tail(components,1);
		all_categories_list[[shortened_category_names[catidx]]]=1;
	}
	#print(shortened_category_names);
	
	# Compute organism totals
	orig_cat_totals=rep(0,num_orig_categories);
	for(cat_idx in 1:num_orig_categories){
		orig_cat_totals[cat_idx]=sum(count_mat[,cat_idx]);
	}
	#print(orig_cat_totals);

	# Compute all totals
	all_total=sum(orig_cat_totals);

	# Normalize cat totals
	normalized=orig_cat_totals/all_total;
	#print(normalized);

	# Store set of cutoffs for each input file
	categories_at_cutoff=list();
	for(cutoff_idx in 1:length(cutoffs)){

		# Apply cutoff 
		cur_cutoff=sprintf("%.2f",cutoffs[cutoff_idx]);
		cat("Applying cutoff: ", cur_cutoff, "\n");
		categories_exceeding_cutoff=logical();
		categories_exceeding_cutoff=normalized>cutoffs[cutoff_idx];

		# Store list of categories for each cutoff
		includeable=shortened_category_names[categories_exceeding_cutoff];
		categories_at_cutoff[[cur_cutoff]]=list();
		for(categories in includeable){
			categories_at_cutoff[[cur_cutoff]][[categories]]=1;
		}

	}

	loaded_data[[InputFileName]]=categories_at_cutoff;
	
	arg_count=arg_count+1;
}

num_inputs=arg_count-1;
cat("Num files: ", num_inputs, "\n");
#print(loaded_data);


venn_counts_list=list();
overlap_types=c("A","B","C","AB","AC","BC","ABC");
for(cutoff in cutoffs_str){
	venn_counts_list[[cutoff]]=list();
	for(overlap_type in overlap_types){
		venn_counts_list[[cutoff]][[overlap_type]]=0;
	}
}

list_keys=names(all_categories_list);
# Compute venn memberships
for(cutoff in cutoffs_str){
	#cat("Working on cutoff: ", cutoff, "\n");
	for(catname in list_keys){
		#cat("  ", catname, "\n");
		presence=numeric(0);
		for(filename in names(loaded_data)){
			presence=c(presence,!is.null(loaded_data[[filename]][[cutoff]][[catname]]));
			cat("     ", filename, " ", presence, "\n");
		}
		if(length(presence)==2){
			presence=c(presence,0);
		}
		#print(presence);

		if(all(presence==c(1,0,0))){
			venn_counts_list[[cutoff]][["A"]]=venn_counts_list[[cutoff]][["A"]]+1;
		}else if(all(presence==c(0,1,0))){
			venn_counts_list[[cutoff]][["B"]]=venn_counts_list[[cutoff]][["B"]]+1;
		}else if(all(presence==c(0,0,1))){
			venn_counts_list[[cutoff]][["C"]]=venn_counts_list[[cutoff]][["C"]]+1;
		}else if(all(presence==c(1,1,0))){
			venn_counts_list[[cutoff]][["AB"]]=venn_counts_list[[cutoff]][["AB"]]+1;
		}else if(all(presence==c(1,0,1))){
			venn_counts_list[[cutoff]][["AC"]]=venn_counts_list[[cutoff]][["AC"]]+1;
		}else if(all(presence==c(0,1,1))){
			venn_counts_list[[cutoff]][["BC"]]=venn_counts_list[[cutoff]][["BC"]]+1;
		}else if(all(presence==c(1,1,1))){
			venn_counts_list[[cutoff]][["ABC"]]=venn_counts_list[[cutoff]][["ABC"]]+1;
		}else if(all(presence==c(0,0,0))){
			# No op.
		}else{
			cat("ERROR! Unknown presence...\n");
			print(presence);
			q(status=-1);
		}
	}	
	#print( venn_counts_list[[cutoff]]);
}

clean_names=sub(".summary_table.xls$","",names(loaded_data));
clean_names=sub(".xls$","",clean_names);

outfname=paste(paste(clean_names, collapse="_"), "venn.pdf", sep=".");

pdf(outfname, height=8.5, width=11);


cat("\n");
for(cutoff in cutoffs_str){
	cat("Graphing Venn Diagram at cutoff = ", cutoff, "\n");
	#print(venn_counts_list[[cutoff]]);
	plot_venn(venn_counts_list[[cutoff]],clean_names, cutoff);
	cat("Done.\n");
}


dev.off();

###############################################################################


writeLines("Done.\n")

print(warnings());
q(status=0)
