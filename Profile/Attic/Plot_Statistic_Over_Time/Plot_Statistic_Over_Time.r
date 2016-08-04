#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');

params=c(
	"treatment_info_file", "t", 1, "character",
	"date_column", "d", 1, "numeric",
	"treatment_column", "c", 1, "numeric",

	"statistics_info_file", "v", 1, "character",
	"statistics_column", "s", 1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <input treatment table (tab separated)>\n",
	"	  -d <column in treatment table that contains the date/time information>\n",
	"	  -c <column in treatment table that contains the control/treatment information>\n",
	"\n",
	"	-v <input statistics table (comma separated)>\n",
	"	  -s <column in statistics table that contains the statistic to compare>\n",
	"\n",
	"Column names count from 0.  So the first non-sample ID column is 1\n",
	"\n",
	"\n");

if(
	!length(opt$treatment_info_file)||
	!length(opt$date_column)||
	!length(opt$treatment_column)||
	!length(opt$statistics_info_file)||
	!length(opt$statistics_column)
){
	cat(usage);
	q(status=-1);
}

TreatmentFileName=opt$treatment_info_file;
DateColumn=opt$date_column;
TreatmentColumn=opt$treatment_column;
StatisticsFileName=opt$statistics_info_file;
StatisticsColumn=opt$statistics_column;


##############################################################################

cat("Treatment File Name: ", TreatmentFileName, "\n");
cat("Statistics File Name: ", StatisticsFileName, "\n");
cat("\n");

##############################################################################

# Load inputs
treatment_matrix=as.matrix(read.table(TreatmentFileName, sep="\t", header=T, row.names=1));
treatment_col_names=colnames(treatment_matrix);
#print(treatment_matrix);

statistics_matrix=as.matrix(read.table(StatisticsFileName, sep=",", header=T, row.names=1));
statistics_col_names=colnames(statistics_matrix);
statistics_sample_names=rownames(statistics_matrix);
#print(statistics_matrix);

##############################################################################

cat("Date/Time column: ", treatment_col_names[DateColumn], "\n");
cat("Treatment column: ", treatment_col_names[TreatmentColumn], "\n");
cat("Statistics column: ", statistics_col_names[StatisticsColumn], "\n");
cat("\n");

OutputFileNameRoot=paste(gsub("\\.csv", "", StatisticsFileName), ".", statistics_col_names[StatisticsColumn], sep="");

##############################################################################
# Extract only needed columns

date_categories=sort(unique(treatment_matrix[,DateColumn]));
date_categories=date_categories[nchar(date_categories)>0];
cat("Date Partitions: ", date_categories, "\n");

treatment_categories=unique(treatment_matrix[,TreatmentColumn]);
treatment_categories=treatment_categories[nchar(treatment_categories)>0];
cat("Treatment types: ", treatment_categories, "\n");

cat("\n");

# Pull out sample names by treatment and time values
treatment_list=length(treatment_categories);
for(trt in treatment_categories){
	cat("Extracting: '", trt, "'\n", sep="");
	treatment_list[[trt]]=list(length(date_categories));
	for(date in date_categories){
		cat("  Extracting: '", date, "'\n", sep=""); 
		date_match=(treatment_matrix[,DateColumn]==date)
		trt_match=(treatment_matrix[,TreatmentColumn]==trt)
		both_match_idx=which(date_match & trt_match);
		#print(both_match_idx);	
		treatment_list[[trt]][[date]]=names(both_match_idx);
	}
}

#print(treatment_list);

##############################################################################


sample_values=list();
stat_min=NA;
stat_max=NA;
for(trt in treatment_categories){
	for(date in date_categories){
		cat("Treatment: ", trt, "Date: ", date, "\n");
		sample_ids=treatment_list[[trt]][[date]];

		intersected=intersect(sample_ids, statistics_sample_names);	
		missing=setdiff(sample_ids, statistics_sample_names);
		if(length(missing)>0){
			cat("\nMissing Sample IDs from Statistics Table: \n");
			print(missing);
			cat("\n");
		}
	
		sample_values[[trt]][[date]]=statistics_matrix[intersected, StatisticsColumn];

		if(is.na(stat_min)){
			stat_min=sample_values[[trt]][[date]];
		}else{
			stat_min=min(stat_min, sample_values[[trt]][[date]]);
		}

		if(is.na(stat_max)){
			stat_max=sample_values[[trt]][[date]];
		}else{
			stat_max=max(stat_max, sample_values[[trt]][[date]]);
		}
		
	}
}

print(sample_values);
print(stat_min);
print(stat_max);
dates=sort(as.numeric(date_categories));
date_min=min(dates);
date_max=max(dates);
colors=c("blue", "dark red");

pdf(paste(OutputFileNameRoot, ".pdf", sep=""), height=8.5, width=11);

plot(0,0, ylim=c(stat_min, stat_max), xlim=c(date_min, date_max), 
	ylab=statistics_col_names[StatisticsColumn], xlab=treatment_col_names[DateColumn],
	type="n");

legend(date_max*.75, stat_max, legend=treatment_categories, fill=colors);


col_idx=1;
for(trt in treatment_categories){
	for(date in date_categories){
		date_val=as.numeric(date);
		stats_array=sample_values[[trt]][[date]];
		num_stats=length(stats_array);
		points(rep(date_val, num_stats), stats_array, col=colors[col_idx], pch=1, cex=.5);
		points(date_val, median(stats_array), col=colors[col_idx], pch=3, cex=2);
	}
	col_idx=col_idx+1;
}

##############################################################################

cat("Done.\n")
dev.off();

q(status=0)
