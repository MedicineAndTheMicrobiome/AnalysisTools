#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('biomformat');

# https://bioconductor.org/packages/release/bioc/html/biomformat.html

params=c(
	"biom_file", "b", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-b <input biom file name>\n",
	"\n",	
	"Converts a biom file into a summary table.\n",
	"\n");

if(!length(opt$biom_file)){
	cat(usage);
	q(status=-1);
}

BiomFilename=opt$biom_file;
OutputFileNameRoot=gsub("\\.biom$", "", BiomFilename);
SummaryTableOut=paste(OutputFileNameRoot, ".summary_table.tsv", sep="");

###############################################################################

cat("\n")
cat("Input Biom File Name: ", BiomFilename, "\n");
cat("Output Summary Table: ", SummaryTableOut, "\n");
cat("\n");

###############################################################################

biom_data=read_biom(BiomFilename);

fields=names(biom_data);

cat("Fields in the biom file:\n");
print(fields);

cat("\n");
for(f in fields){
	cat("Field: ", f, "\n");
	if(any(f==c("data", "rows", "columns"))){
		cat("[Too much to show.]\n");
	}else{
		print(biom_data[[f]]);
	}
	cat("\n");
}


num_categories=biom_data$shape[1];
num_samples=biom_data$shape[2];

cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

otu_id_width=ceiling(log10(num_categories));

cat("Extracting categories...\n");

category_names=character(num_categories);

for(c in 1:num_categories){
	cat_info=biom_data$rows[c];
	classifications=cat_info[[1]]$metadata$taxonomy;

	otu_id=sprintf(paste("_otu%0", otu_id_width, "i", sep=""), as.numeric(cat_info[[1]]$id));

	taxa_str=paste(classifications, collapse=";");

	for(t in 1:7){
		txt=classifications[t];
		txt=gsub("^.__", "", txt);

		if(txt=="" || txt=="__"){
			txt=classifications[t-1];
		}else{
			txt=gsub("\\[", "", txt);
			txt=gsub("\\]", "", txt);
			txt=gsub("-", "_", txt);

		}

		if(t==7){
			classifications[t]=paste(txt, otu_id, sep="");
		}else{
			classifications[t]=txt;;
		}
	}

	category_names[c]=paste(classifications, collapse=";");
}

cat("Extracting samples...\n");

sample_names=character(num_samples);
for(r in 1:num_samples){
	samp_info=biom_data$columns[r];
	sample_names[r]=samp_info[[1]]$id;
}


cat("Extracting data...\n");

counts=matrix(NA, nrow=num_samples, ncol=num_categories);
rownames(counts)=sample_names;
colnames(counts)=category_names;

for(c in 1:num_categories){
	cat_data=(biom_data$data[[c]]);
	s_name=names(cat_data);
	counts[s_name, c]=cat_data;
}


# Output
cat("\nWriting New Matrix...\n");

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

write_summary_file(counts, SummaryTableOut);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
