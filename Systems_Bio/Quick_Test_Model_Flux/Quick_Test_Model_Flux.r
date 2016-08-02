#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
		"model", "m", 1, "character",
		"flux", "f", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <model name>\n",
		"	[-f <output flux filename>]\n",
		"\n",
		"This script will read in the model, and compute the \n",
		"flux through the model.\n",
		"\n",
		"The model name: Thiobacillus_Denitrificans, for example,\n",
		"would read in the model with the following 3 components:\n",
		"\tThiobacillus_Denitrificans_react.tsv\n",
		"\tThiobacillus_Denitrificans_met.tsv\n",
		"\tThiobacillus_Denitrificans_desc.tsv\n",
		"\n");

if(!length(opt$model) ){
	cat(usage);
	q(status=-1);
}

library(sybil);
library(glpkAPI);
cat("\n\n");

FluxFname="";
if(length(opt$flux)){
	FluxFname=opt$flux;
	cat("Output Flux Filename: ", FluxFname, "\n", sep="");
}

Model=opt$model;

# Clean the name just in case more than the prefix was specified.
Model=gsub("_[a-z]+\\.tsv$","",Model);
cat("Model: ", Model, "\n", sep="");

###############################################################################

test_model=function(model){
	opt_res=optimizeProb(model, algorithm="fba", retOptSol=TRUE);

	objective_flux=attributes(opt_res)$lp_obj;
	reaction_flux=(attributes(attributes(opt_res)$fluxdist)$fluxes);

	optim_out=list();
	optim_out$objective_flux=objective_flux;
	optim_out$reaction_flux=reaction_flux;
	
	return(optim_out);
}

###############################################################################

# Load model
path_arr=strsplit(Model, "/")[[1]];

model_prefix=tail(path_arr, 1);
model_path=paste(head(path_arr, length(path_arr)-1), collapse="/");

if(model_path==""){
	model_path=paste(".",model_path, sep="");
}

met_info=file.info(paste(Model, "_met.tsv", sep=""));
react_info=file.info(paste(Model, "_react.tsv", sep=""));
desc_info=file.info(paste(Model, "_desc.tsv", sep=""));

cat("\n");
cat("Last Modified:\n");
cat("  Description: ", format(as.POSIXct(desc_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"), "\n");
cat("  Metabolites: ", format(as.POSIXct(met_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"),  "\n");
cat("    Reactions: ", format(as.POSIXct(react_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"), "\n");
cat("\n");

cat("  Model Path: ", model_path, "\n");
cat("Model Prefix: ", model_prefix, "\n");
cat("\n");


model=readTSVmod(fpath=model_path, prefix=model_prefix, suffix="tsv");

#print(met_comp(model));
#print(S(model));

###############################################################################

get_met_compartments=function(model){
	avail_compartments=mod_compart(model);
	if(length(avail_compartments)==1){
		avail_compartments=strsplit(avail_compartments, ",")[[1]];
	}
	num_avail_compart=length(avail_compartments);
	met_ids=met_id(model);	
	
	compart_ids=character();
	for(i in 1:num_avail_compart){
		ix=grep(paste("\\[", avail_compartments[i], "\\]", sep=""), met_ids);
		compart_ids[ix]=avail_compartments[i];
	}
	
	return(compart_ids);
}

get_rxn_compartments=function(model){
	
	met_compartments=get_met_compartments(model);

	S_mat=S(model);
	num_rxns=ncol(S_mat);

	rxn_comp=character();

	for(i in 1:num_rxns){
		met=S_mat[,i];
		met_ix=which(met!=0);
		met_comp=met_compartments[met_ix];
		uniq_met_comp=sort(unique(met_comp));
		rxn_comp[i]=paste(uniq_met_comp, collapse="/");
	}
	
	return(rxn_comp);	

}

build_rxn_string=function(model){
	S_mat=S(model);
	num_rxns=ncol(S_mat);
	met_ids=met_id(model);
	
	rxn_str=list();
	for(i in 1:num_rxns){
		lhs=which(S_mat[,i]<0);
		rhs=which(S_mat[,i]>0);

		# Generate coefficient string
		lhs_coef_str=paste("(", abs(S_mat[lhs, i]), ")", sep="");
		rhs_coef_str=paste("(", abs(S_mat[rhs, i]), ")", sep="");

		# Remove coefficient string if it's just (1)
		lhs_coef_str=gsub("\\(1\\)", "", lhs_coef_str);
		rhs_coef_str=gsub("\\(1\\)", "", rhs_coef_str);

		# Append coefficient string to metabolite ID
		lhs_met=paste(lhs_coef_str, met_ids[lhs]);
		rhs_met=paste(rhs_coef_str, met_ids[rhs]);

		# Concatenate metabolites IDs together
		lhs_str=paste(lhs_met, collapse=" + ");
		rhs_str=paste(rhs_met, collapse=" + ");

		# Remove leading and double spaces
		lhs_str=gsub("  +", " ", lhs_str);
		rhs_str=gsub("  +", " ", rhs_str);
		lhs_str=gsub("^ +", "", lhs_str);
		rhs_str=gsub("^ +", "", rhs_str);

		rxn_str[[i]]=list();
		rxn_str[[i]][["LHS"]]=lhs_str;
		rxn_str[[i]][["RHS"]]=rhs_str;
	}

	return(rxn_str);

}

###############################################################################

optim_out=test_model(model);

MIN_TOL=0.0001;

succ=logical();
if(optim_out$objective_flux>MIN_TOL){
	succ=T;
}else{
	succ=F;
}

if(succ){
	cat("\nModel produced flux:  ", optim_out$objective_flux, "\n");
}else{
	cat("\nModel did not produce flux: ", optim_out$objective_flux, "\n");
}

rxn_id=react_id(model);
rxn_names=react_name(model);
rxn_comp=get_rxn_compartments(model);
rxn_strs=build_rxn_string(model);

if(succ && FluxFname!=""){
	fh=file(paste(FluxFname, ".tsv", sep=""), "w");
	cat(file=fh, paste(c("Id","Description","Flux","Compartment","Equation"), collapse="\t"), "\n", sep="");

	num_rxns=nrow(optim_out$reaction_flux);
	for(i in 1:num_rxns){

		flux=optim_out$reaction_flux[i,1];

		if(abs(flux)>1e-10){
			if(flux>0){
				dir="-->";
			}else if(flux<0){
				dir="<--";
			}else{
				dir="-|-";
			}


			rxn_eq=paste(rxn_strs[[i]][["LHS"]], dir, rxn_strs[[i]][["RHS"]], sep=" ");

			if(flux==0){
				rxn_eq="";
			}		

			cat(file=fh, 
				rxn_id[i], rxn_names[i], flux, rxn_comp[i], 
				rxn_eq, sep="\t");
			cat(file=fh, "\n");
		}
	}
}

###############################################################################

print(warnings());
cat("Done.\n");


























