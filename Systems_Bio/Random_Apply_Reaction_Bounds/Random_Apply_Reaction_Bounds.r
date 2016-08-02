#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
		"model", "m", 1, "character",
		"reaction", "r", 1, "character",
		"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <model name>\n",
		"	-r <reaction list>\n",
		"	-o <output root>\n",
		"\n",
		"This script will read in the model, and then go through\n",
		"the reaction list, sequentially apply and test each\n",
		"reaction in the reaction list.\n",
		"\n",
		"The model name: Thiobacillus_Denitrificans, for example,\n",
		"would read in the model with the following 3 components:\n",
		"\tThiobacillus_Denitrificans_react.tsv\n",
		"\tThiobacillus_Denitrificans_met.tsv\n",
		"\tThiobacillus_Denitrificans_desc.tsv\n",
		"\n");

if(!length(opt$model) || 
	!length(opt$reaction) ||
	!length(opt$output)
){
	cat(usage);
	q(status=-1);
}

library(sybil);
library(glpkAPI);

ReactionList=opt$reaction;
Model=opt$model;
Output=opt$output;

cat("Model: ", Model, "\n");
cat("Reaction List: ", ReactionList, "\n");

###############################################################################

load_list=function(fname){
	rxn_attr=read.delim(file=fname, sep="\t", header=T, row.names=1, check.names=F);
	#print(rxn_attr);
	return(rxn_attr);
}

###############################################################################

change_bounds=function(model, idx, lb, ub){
	lbs=lowbnd(model);
	ubs=uppbnd(model);
	lbs[idx]=lb;
	ubs[idx]=ub;
	lowbnd(model)=lbs;
	uppbnd(model)=ubs;
	return(model);
}

apply_reaction=function(model, rxn_idx, direction, reversible){

	# direction: 1 forward, -1 reverse
	# reversible: 0 not reversible, 1 reversible

	if(direction==-1){
		lb=-1000;
		ub=0;
	}else{
		lb=0;
		ub=1000;
	}

	if(reversible==1){
		lb=-1000;
		ub=1000;
	}

	new_model=change_bounds(model, rxn_idx, lb, ub);	

	return(new_model);
}

test_model=function(model){
	MIN_TOL=0.0001;

	opt_res=optimizeProb(model, algorithm="fba", retOptSol=TRUE);
	flux=attributes(opt_res)$lp_obj;
	
	succ=logical();
	if(flux>MIN_TOL){
		succ=T;
	}else{
		succ=F;
	}
	return(succ);
}


###############################################################################

# Load reaction list
rxn_attr=load_list(ReactionList);
num_reactions=nrow(rxn_attr);
cat("Num reactions to test: ", num_reactions, "\n");

# Load model
path_arr=strsplit(Model, "/")[[1]];
model_prefix=tail(path_arr, 1);
model_path=paste(head(path_arr, length(path_arr)-1), collapse="/");

cat("Model Path: ", model_path, "\n");
cat("Model Prefix: ", model_prefix, "\n");
model=readTSVmod(fpath=model_path, prefix=model_prefix, suffix="tsv");

###############################################################################

num_trials=30;
#num_trials=num_reactions;

react_ids=react_id(model)
test_react_names=rownames(rxn_attr);

killer_rxns=character(num_trials);

for(t_ix in 1:num_trials){

	modified_model=model;
	sequence=sample(1:num_reactions, num_reactions, replace=F);
	cat("[", t_ix, "] Starting new sequence.\n", sep="");
	
	fully_succ=T;
	for(i in sequence){
		rxn_info=rxn_attr[i,];
		cur_react_id=test_react_names[i];

		cat("  Modifying Reaction ID: ", cur_react_id, "\n", sep="");
		rxn_idx=which(react_ids==cur_react_id);
		
		modified_model=apply_reaction(modified_model, rxn_idx, direction=rxn_info[2], reversible=rxn_info[3]);
		succ=test_model(modified_model);
		if(!succ){
			killer_rxns[t_ix]=cur_react_id;
			fully_succ=F;		
			break;
		}
	}

	if(fully_succ){
		cat("Fully Successful.\n");
		break;
	}	
}

###############################################################################

# Output counts/reaction for reactions that killed the flux
killer_table=table(killer_rxns);
killer_names=names(killer_table);
killer_counts=as.vector(killer_table);

fh=file(paste(Output, ".killer_rxns.txt", sep=""), "w");

num_killer_rxn=length(killer_counts);
for(i in 1:num_killer_rxn){
	cat(file=fh, killer_names[i], "\t", killer_counts[i], "\n", sep="");
}

close(fh);

###############################################################################


print(warnings());
cat("Done.\n");


























