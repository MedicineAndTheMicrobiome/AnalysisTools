#!/usr/bin/Rscript

test_code=T;

sample_sheet_validator=function(df){

	cat("Sample Sheet Validator Started...\n");
	BC_PLATEMAP="barcode_plateloc.tsv";

	first_col=df[,1];	

	# get experiment name:
	cat("Getting Experiment Name...\n");
	exp_name_ix=which(first_col=="Experiment Name");
	experiment_name=df[exp_name_ix, 2];

	# get experiment date:
	cat("Getting Experiment Date...\n");
	exp_date_ix=which(first_col=="Date");
	date=df[exp_date_ix, 2];

	# get investigator name:
	cat("Getting Investigator Name...\n");
	
	investigator_name_ix=(first_col=="Investigator Name");
	if(all(!investigator_name_ix)){
		investigator_name="<Not Specified>";
	}else{
		investigator_name=df[investigator_name_ix, 2];
	}

	# get [Data] and Sample_ID position
	cat("Getting Data Line Number...\n");
	data_ix=which(first_col=="[Data]");
	sample_id_header_ix=which(first_col=="Sample_ID");

	if(data_ix==(sample_id_header_ix-1)){
		sample_begin_ix=sample_id_header_ix+1;
	}else{
		return("Error:  Could not find [Data] and Sample_ID row.\n");
	}

	cat("Sample information detected to start on line: ", sample_begin_ix, "\n"); 

	df_rows=nrow(df);
	df_cols=ncol(df);
	
	#----------------------------------------------------------------------
	# Extract sample information:
	samp_mat=as.matrix(df[sample_begin_ix:df_rows,,drop=F]);
	hdr=as.matrix(df[sample_id_header_ix,]);
	colnames(samp_mat)=hdr;

	# Confirm column names are as expected:
	expected_columns=c("Sample_ID", "Sample_Name", "index", "I7_Index_ID")
	missing_columns=setdiff(expected_columns,hdr);
	if(length(missing_columns)>0){
		miss_col_list=paste(missing_columns, collapse=", ");
		err=paste("Error: Could not find columns: ", miss_col_list, "\n", sep="");
		return(err);
	}
	

	#print(samp_mat)
	num_samples=nrow(samp_mat);
	samp_ids=samp_mat[,"Sample_ID"];
	samp_names=samp_mat[,"Sample_Name"];
	barcodes=samp_mat[,"index"];
	plate_map=samp_mat[,"I7_Index_ID"];

	cat("Starting individual checks...\n");

	#----------------------------------------------------------------------
	# Look for dup barcodes
	cat("Looking for duplicate barcodes...\n");
	bc_uniq=unique(barcodes);
	num_uniq_bc=length(bc_uniq);
	bc_are_uniq=T;
	if(num_uniq_bc>num_samples){
		bc_uniq_msg="ERROR: More barcodes than samples";
	}else if(num_uniq_bc<num_samples){
		counts=table(barcodes);
		dup_ids=paste(names(counts[counts>1]), collapse=",");
		bc_uniq_msg=c("ERROR: Barcodes are not unique:", dup_ids);
	}else{
		bc_uniq_msg="OK: Barcodes look unique";
	}

	#----------------------------------------------------------------------
	# Look for invalid barcodes characters;
	cat("Looking for invalid barcode characters...\n");
	inval_bc_msg=c();
	for(b in barcodes){
		out=gsub("[ATGC]", "", b);
		if(nchar(out)>0){
			inval_bc_msg=c(inval_bc_msg, b);	
		}
	}
	if(length(inval_bc_msg)>0){
		inval_bc_msg=c("ERROR: Invalid barcodes detected:", inval_bc_msg);
	}else{
		inval_bc_msg=c("OK: Barcode sequences look good");
	}

	#----------------------------------------------------------------------
	# Look for inconsistent bc lengths
	cat("Looking for inconsistent barcode lengths...\n");
	incons_bc_len_msg=c();
	bc_lengths=numeric();
	for(i in 1:num_samples){
		bc_lengths[i]=nchar(barcodes[i]);
	}
	uniq_lengths=unique(bc_lengths);
	if(length(uniq_lengths)>1){
		incons_bc_len_msg=c("ERROR: Barcode lengths not consistent:", 
			paste(uniq_lengths, collapse=", "));
	}else{
		incons_bc_len_msg=c(paste("OK: All barcodes the same length: ", uniq_lengths));
	}

	#----------------------------------------------------------------------
	# Look for duplicate sample IDs
	cat("Looking for duplicate sample IDs...\n");
	uniq_samp_ids=unique(samp_names);
	num_uniq_samp_ids=length(uniq_samp_ids);
	uniq_samp_ids_msg=character();
	if(num_uniq_samp_ids<num_samples){
		counts=table(samp_names);
		dup_ids=paste(names(counts[counts>1]), collapse=",");
		uniq_samp_ids_msg=c("ERROR: Sample IDs are not unique:", dup_ids);
	}else{
		uniq_samp_ids_msg="OK: All Sample IDs are unique";
	}

	#----------------------------------------------------------------------
	# Confirm sample ids match sample names
	cat("Confirming that Sample IDs match Sample Names...\n");
	pairs_match=(samp_ids==samp_names);
	if(all(pairs_match)){
		samp_id_and_names_match_msg="OK: All Sample IDs match the Sample Names";
	}else{
		samp_id_and_names_match_msg=c("ERROR: Samples IDs do not match Sample Names",
			paste(samp_ids[!pairs_match],"/",samp_names[!pairs_match]));
	}

	#----------------------------------------------------------------------
	# Check sample IDs for lowercase characters
	cat("Checking Sample IDs for lowercase characters...\n");
	uc_samp_ids=toupper(samp_ids);
	non_uc_ix=!(uc_samp_ids==samp_ids)
	if(any(non_uc_ix)){
		mixcase_msg=c("ERROR: Sample ID's are not all in uppercase.",
			paste(samp_ids[non_uc_ix], collapse=", ")
		);
	}else{
		mixcase_msg=c("OK: All Sample IDs have consistent uppercase.");
	}

	
	#----------------------------------------------------------------------
	# Check sample IDs for invalid characters
	cat("Checking Sample IDs for invalid characters...\n");
	bad_samp_ids=character();
	bad_char=character();
	for(i in 1:num_samples){
		inval_char=gsub("[A-Z0-9_]", "", samp_ids[i]);
		if(nchar(inval_char)>0){
			bad_samp_ids=c(bad_samp_ids, samp_ids[i]);
			bad_char=c(bad_char, inval_char);
		}
	}
	if(length(bad_samp_ids)>0){
		bad_char_msg=c("ERROR: Bad characters found in Sample IDs",
			paste(bad_samp_ids, " (", bad_char, ")", sep="")
		);
	}else{
		bad_char_msg="OK: All Sample IDs have valid characters.";
	}


	#----------------------------------------------------------------------
	cat("Checking for underscore at end of sample ID...\n");
	bad_samp_ids=character();
	for(i in 1:num_samples){
		if(length(grep("_$", samp_ids[i]))>0){
			bad_samp_ids=c(bad_samp_ids, samp_ids[i]);
		}
	}
	if(length(bad_samp_ids)>0){
		bad_char_endunderscore_msg=c("ERROR: Underscore found at end of Sample IDs",
			paste(bad_samp_ids, collapse=", ")
		);
	}else{
		bad_char_endunderscore_msg="OK: All Sample IDs don't have underscores at end";
	}

	#----------------------------------------------------------------------
	# Report on project codes	
	cat("Reporting on Project Codes...\n");
	project_codes=character(num_samples);

	splits=strsplit(samp_ids, "_");
	for(i in 1:num_samples){
		project_codes[i]=splits[[i]][1];	
	}
	prj_tab=table(project_codes);
	prj_mat=matrix(prj_tab, nrow=length(prj_tab), ncol=1);
	rownames(prj_mat)=names(prj_tab);
	colnames(prj_mat)="counts";
	pcode_tab=capture.output(print(prj_mat, quote=F));
	pcode_tot=sum(prj_tab);
	pcode_tab=c(pcode_tab, 
		"--------------------", 
		paste("Total", pcode_tot));

	#----------------------------------------------------------------------
	# project code lengths	
	cat("Checking Project Code lengths...\n");
	uniq_pcodes=names(prj_tab);
	non_stnd_len_pcodes=c();
	for(i in 1:length(uniq_pcodes)){
		if(nchar(uniq_pcodes[i])!=4){
			non_stnd_len_pcodes=c(non_stnd_len_pcodes, uniq_pcodes[i]);	
		}
	}
	if(length(non_stnd_len_pcodes)>0){
		non_stnd_len_pcodes_msg=
			c("ERROR: Some projects codes are not 4 characters.",
			non_stnd_len_pcodes);
	}else{
		non_stnd_len_pcodes_msg="OK: Project codes have a consistent length";
	}

	#----------------------------------------------------------------------
	# Check for XX's placeholders
	cat("Checking for more than 3 consecutive X's in Sample ID...\n");
	bad_samp_ids=character();
	
	for(i in 1:num_samples){
		cur_id=samp_ids[i];
		if(length(grep("XXX", cur_id))){
			bad_samp_ids=c(bad_samp_ids, cur_id);
		}
	}
	if(length(bad_samp_ids)>0){
		residual_xxx_placeholder_msg=
			c("ERROR: Sample IDs with XXX's identified.",
				bad_samp_ids);
	}else{
		residual_xxx_placeholder_msg="OK: No Sample IDs with XXX's identified";
	}


	#----------------------------------------------------------------------
	# cross check Sample_Project (which is the plate_id/well for the barcode) with Barcode

	# Load barcodes
	plate_map_ref=as.matrix(read.table(BC_PLATEMAP, header=T, sep=","));	
	num_pm_entries=nrow(plate_map_ref);
	pm_list=list();
	for(i in 1:num_pm_entries){
		key=plate_map_ref[i,1];
		if(key==""){
			next;
		}
		pm_list[[key]]=plate_map_ref[i,2];
	}
	cat("Number of Plate Map Entries Read from ", BC_PLATEMAP, ": ", num_pm_entries, "\n", sep="");
	#print(pm_list);

	cat("Checking Plate Map consistency...\n");
	no_info_cnt=0;
	no_info_lst=c();
	mism_cnt=0;
	mism_info_lst=c();
	for(i in 1:num_samples){
		pm=plate_map[i];
		bc=barcodes[i];

		exp_bc=pm_list[[pm]];

		if(length(exp_bc)==0){
			no_info_lst=c(no_info_lst,
				paste("Line: ", i+(sample_begin_ix-1), 
					", has no plate map information: ", pm, " / ", bc, sep=""));
			no_info_cnt=no_info_cnt+1;
		}else if(exp_bc!=bc){
			mism_info_lst=c(mism_info_lst,
			paste("Line: ", i+(sample_begin_ix-1), ", PlateMap/Barcode Mismatch:", pm, "/", bc, sep=""),
			paste(" Expected: ", pm, "/", pm_list[[pm]], sep=""));
			mism_cnt=mism_cnt+1;
		}
	}

	if(no_info_cnt>0){
		plate_info_msg=c("ERROR: Missing plate info:", no_info_lst);
	}else{
		plate_info_msg="OK: Plate Info and Barcodes Accounted For ";
	}

	if(mism_cnt>0){
		plate_match_msg=c("ERROR: Mismatching plate info:", mism_info_lst);
	}else{
		plate_match_msg="OK: Plate Info and Barcodes Match";
	}

	#----------------------------------------------------------------------
	cat("Generating output lines...\n");	
	msg_arr=c(	
		paste("Experiment Name: ", experiment_name, sep=""),
		paste("Investigator Name: ", investigator_name, sep=""),
		paste("Date: ", date, sep=""),
		"",
		paste("Overall Sheet Dimensions: ", df_rows, " r x ", df_cols, " c", sep=""),
		paste("Number of Samples Detected: ", num_samples, sep=""),
		"",
		bc_uniq_msg,
		inval_bc_msg,
		incons_bc_len_msg,
		uniq_samp_ids_msg,
		samp_id_and_names_match_msg,
		mixcase_msg,
		bad_char_msg,
		bad_char_endunderscore_msg,
		residual_xxx_placeholder_msg,
		non_stnd_len_pcodes_msg,
		plate_info_msg,
		plate_match_msg,
		"",
		"Project Summary:",
		pcode_tab,
		""
	);	

	cat("Exiting Validator...\n");
	return(msg_arr);
}


if(test_code && F){

	cat("-----------------------------------------------------------------\n");

	df=read.csv("example/20200819_a_saliva_samplesheet.csv", header=F);
	msgs=sample_sheet_validator(df);
	cat(paste(msgs, collapse="\n"));
}


if(test_code && F){
	cat("-----------------------------------------------------------------\n");
	df=read.csv("example/0159_20200131_samplesheet.wErrors.csv", header=F);
	msgs=sample_sheet_validator(df);
	cat(paste(msgs, collapse="\n"));

	cat("-----------------------------------------------------------------\n");
}

print(warnings());


