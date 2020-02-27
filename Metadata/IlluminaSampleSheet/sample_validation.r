#!/usr/bin/Rscript

# Look for duplicate barcodes
# Invalid barcodes
# Duplicate sample IDs
# Sample IDs with bad names
# Missing 4 digit project code



sample_sheet_validator=function(df){

	first_col=df[,1];	

	# get experiment name:
	exp_name_ix=which(first_col=="Experiment Name");
	experiment_name=df[exp_name_ix, 2];

	# get experiment date:
	exp_date_ix=which(first_col=="Date");
	date=df[exp_date_ix, 2];

	# get [Data] and Sample_ID position
	data_ix=which(first_col=="[Data]");
	sample_id_header_ix=which(first_col=="Sample_ID");

	if(data_ix==(sample_id_header_ix-1)){
		sample_begin_ix=sample_id_header_ix+1;
	}else{
		return("Error:  Could not find [Data] and Sample_ID row.\n");
	}

	cat("Sample information should start on line: ", sample_begin_ix, "\n"); 

	df_rows=nrow(df);
	df_cols=ncol(df);
	
	# Extract sample information:
	samp_mat=as.matrix(df[sample_begin_ix:df_rows,,drop=F]);
	hdr=as.matrix(df[sample_id_header_ix,]);
	colnames(samp_mat)=hdr;

	print(samp_mat);

	num_samples=nrow(samp_mat);

	# Look for dup barcodes
	bc=samp_mat[,"index"];
	bc_uniq=unique(bc);
	num_uniq_bc=length(bc_uniq);
	bc_are_uniq=T;
	if(num_uniq_bc>num_samples){
		bc_uniq_msg="ERROR: More barcodes than samples";
	}else if(num_uniq_bc<num_samples){
		counts=table(bc);
		dup_ids=paste(names(counts[counts>1]), collapse=",");
		bc_uniq_msg=c("ERROR: Barcodes are not unique:", dup_ids);
	}else{
		bc_uniq_msg="OK: Barcodes look unique";
	}

	# Look for invalid barcodes characters;
	inval_bc_msg=c();
	for(b in bc){
		out=gsub("[ATGC]", "", b);
		if(nchar(out)>0){
			inval_bc_msg=c(inval_bc_msg, b);	
		}
	}
	if(length(inval_bc_msg)>0){
		inval_bc_msg=c("ERROR: Invalid barcodes detected", inval_bc_msg);
	}else{
		inval_bc_msg=c("OK: Barcode sequences look good");
	}

	# Look for inconsistent bc lengths
	incons_bc_len_msg=c();
	bc_lengths=numeric();
	for(i in 1:num_samples){
		bc_lengths[i]=nchar(bc[i]);
	}
	uniq_lengths=unique(bc_lengths);
	if(length(uniq_lengths)>1){
		incons_bc_len_msg=c("ERROR: Barcode lengths not consistent:", 
			paste(uniq_lengths, collapse=", "));
	}else{
		incons_bc_len_msg=c(paste("OK: All barcodes the same length: ", uniq_lengths));o
	}


	
	msg_arr=c(	
		paste("Experiment Name: ", experiment_name, sep=""),
		paste("Date: ", date, sep=""),
		"",
		paste("Overall Sheet Dimensions: ", df_rows, " r x ", df_cols, " c", sep=""),
		paste("Number of Samples Detected: ", num_samples, sep=""),
		"",
		bc_uniq_msg,
		inval_bc_msg,
		incons_bc_len_msg,
		""
	);	

}


df=read.csv("example/0159_20200131_samplesheet.csv", header=F);

#print(df);
msgs=sample_sheet_validator(df);
cat(paste(msgs, collapse="\n"));
