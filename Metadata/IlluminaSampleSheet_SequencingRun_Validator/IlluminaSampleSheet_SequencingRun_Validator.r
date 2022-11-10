#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"sequencing_run_directory", "d", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-d <sequencing run directory>\n",
	"\n",
	"This script will look in the sequencing run directory\n",
	"for the sample sheet and then try to find the fastq.gz\n",	
	"file in the Runs directory.\n",
	"\n",
	"Make sure the sample sheet can be found in:\n",
	"<sequencing run directory>/Samplesheet/<>.csv\n",
	"\n",
	"and the sequence data is in:\n",
	"<sequencing run directory>/Run\n",
	"\n");

if(!length(opt$sequencing_run_directory)){
	cat(usage);
	q(status=-1);
}

SeqRunDir=opt$sequencing_run_directory;

cat("\n");
cat("Starting Sample Sheet / Sequencing Run Cross Validator...\n");
cat("\n");
cat("Sequencing Run Directory: ", SeqRunDir, "\n", sep="");
cat("\n");

##############################################################################

get_sample_ids_from_samplesheet=function(sample_sheet_fname){


	#Finding line with [Data] header

	ss_data=readLines(sample_sheet_fname)
	data_row=grepl("^\\[Data\\]", ss_data)
	data_row=which(data_row==TRUE)

	#Basic sanity checks

        if(length(data_row)==0){
                cat("Error: Could not find [Data] in sample sheet.\n");
                quit(status=-1);
        }


	if(!grep("Sample_ID", ss_data[data_row+1])){
               cat("Error: Could not find Sample_ID in sample sheet.\n");
               quit(status=-1);
       }


	#Extracting sampleIDs from sample sheet.

	sample_ids=as.character(ss_data[(data_row+2):NROW(ss_data)]);
	sample_ids=sapply(strsplit(sample_ids, ","), getElement, 1)
	num_samp_ids=length(sample_ids);


	cat("Number of Sample IDs Found in sample sheet: ", num_samp_ids, "\n");

	return(sort(sample_ids));

}

#------------------------------------------------------------------------------

get_sample_ids_from_fastqgz_dir=function(run_dir){

	files=list.files(run_dir, recursive=T);

	# Find fastq.gz files
	targets_ix=grep("\\.fastq\\.gz", files);
	fastq_paths_list=files[targets_ix];

	# Strip path	
	num_fq_files=length(fastq_paths_list);

	fastqgz_list=character();
	for(i in 1:num_fq_files){
		fastqgz_list[i]=tail(strsplit(fastq_paths_list[i], "/")[[1]], 1);
	}

	r1_ix=grep("_S\\d+\\_L001\\_R1_001\\.fastq\\.gz", fastqgz_list);
	r2_ix=grep("_S\\d+\\_L001\\_R2_001\\.fastq\\.gz", fastqgz_list);

	r1_names=fastqgz_list[r1_ix];
	r2_names=fastqgz_list[r2_ix];

	r1_samp_ids=gsub("_S\\d+\\_L001\\_R1_001\\.fastq\\.gz", "", r1_names);
	r2_samp_ids=gsub("_S\\d+\\_L001\\_R2_001\\.fastq\\.gz", "", r2_names);

	#print(r1_samp_ids);
	#print(r2_samp_ids);

	excl_r1=setdiff(r1_samp_ids, r2_samp_ids);
	excl_r2=setdiff(r2_samp_ids, r1_samp_ids);

	if(length(excl_r1)){
		cat("Warning: Missing R2 files for these R1 files:\n");
		print(excl_r1);
	}

	if(length(excl_r2)){
		cat("Warning: Missing R1 files for these R2 files:\n");
		print(excl_r2);
	}

	complete=intersect(r1_samp_ids, r2_samp_ids);
	incomplete=setdiff(union(r1_samp_ids, r2_samp_ids), complete);

	out=list();
	out[["R1_missing"]]=sort(excl_r2);	
	out[["R2_missing"]]=sort(excl_r1);	
	out[["incomplete"]]=sort(incomplete);
	out[["paired"]]=sort(complete);

	return(out);
}

##############################################################################

plot_text=function(strings){
        par(family="Courier");
        par(oma=rep(.1,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }
}

##############################################################################

# Find sample sheet
samplesheet_dir=paste(SeqRunDir, "/Samplesheet", sep="");

cat("Samplesheet Directory: ", samplesheet_dir, "\n");

Samplesheet_dir_contents=list.files(samplesheet_dir, pattern="*\\.csv$");

if(length(Samplesheet_dir_contents)==0){
	cat("Error: Could not find sample sheet in: ", samplesheet_dir, "\n");
	quit(status=-1);
}

if(length(Samplesheet_dir_contents)>1){
	cat("Error: More than one file found that could be sample sheet.\n");
	print(Samplesheet_dir_contents);
	quit(status=-1);
}

cat("\n");
cat("Found Samplesheet: ", Samplesheet_dir_contents, "\n");
cat("\n");

# Get IDs from sample sheet
sample_sheet_path=paste(samplesheet_dir, "/", Samplesheet_dir_contents, sep="");
ss_ids=get_sample_ids_from_samplesheet(sample_sheet_path);

# Get IDs from run directory
fastq_directory=paste(SeqRunDir, "/Run", sep="");
run_ids_rec=get_sample_ids_from_fastqgz_dir(fastq_directory);
run_ids_rec[["paired"]]=sapply(run_ids_rec[["paired"]], function(x){ gsub("-", ".", x);});
run_ids_rec[["R1_missing"]]=sapply(run_ids_rec[["R1_missing"]], function(x){ gsub("-", ".", x);});
run_ids_rec[["R2_missing"]]=sapply(run_ids_rec[["R2_missing"]], function(x){ gsub("-", ".", x);});

# Convert -'s and _'s to .'s to make them consistent
ss_ids_period=sapply(ss_ids, function(x){ gsub("_", ".", x);});
num_ss_ids=length(ss_ids_period);
run_ids_period=run_ids_rec[["paired"]];
num_run_ids=length(run_ids_period);

names(ss_ids_period)=c();
names(run_ids_period)=c();

cat("Samples IDs from Sample Sheet:\n");
print(ss_ids_period);

cat("\n");
cat("Sample IDs with R1 & R2 fastq.gz files:\n");
print(run_ids_period);

missing_fq_files=setdiff(ss_ids_period, run_ids_period);
missing_ss_files=setdiff(run_ids_period, ss_ids_period);

num_missing_fq_files=length(missing_fq_files);
print(missing_ss_files);
num_excess_ss_entries=length(missing_ss_files);

###############################################################################
# Output to PDF file

output_dir=paste(SeqRunDir, "/smpsht_fastq_xval.pdf", sep="");
pdf(output_dir, height=8.5, width=11);

plot_text(c(
	"Samplesheet / Run Fastq File Cross Validator",
	"",
	"",
	paste("Sequencing Run Directory: ", SeqRunDir, sep=""),
	paste("Date/Time Executed: ", date()),
	"",
	"Samplesheet found at: ", 
	paste("  ", sample_sheet_path),
	"",
	"Fastq directory at: ", 
	paste("  ", fastq_directory),
	"",
	"",
	paste("Num Sample IDs in samplesheet: ", num_ss_ids, sep=""),
	paste("Num Sample IDs by paired fastq files: ", num_run_ids, sep=""),
	"",
	"",
	paste("Num Missing Fastq Samples IDs: ", num_missing_fq_files, sep=""),
	paste("Num Excess Sample Sheet IDs: ", num_excess_ss_entries, sep=""),
	"",
	"",
	{if(num_ss_ids!=num_run_ids){"Error: Sample ID count mismatch."}}
	
));


plot_text(c(
	"Sample IDs missing fastq.gz files:",
	"",
	{if(length(missing_fq_files)){missing_fq_files}else{ "(none)"}}
));

plot_text(c(
	"Sample IDs missing from sample sheet:",
	"",
	{if(length(missing_ss_files)){missing_ss_files}else{"(none)"}}
));

plot_text(c(
	"Missing R1 Fastq Files:",
	"",
	{if(length(run_ids_rec[["R1_missing"]])){run_ids_rec[["R1_missing"]]}else{"(none)"}}
));

plot_text(c(
	"Missing R2 Fastq Files:",
	"",
	{if(length(run_ids_rec[["R2_missing"]])){run_ids_rec[["R2_missing"]]}else{"(none)"}}
));

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);

