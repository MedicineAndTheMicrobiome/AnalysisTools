#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
		"model",    "m", 1, "character",
		"def_flux", "d", 1, "character",
		"ex_rxn_A", "A", 1, "character",
		"ex_rxn_B", "B", 1, "character",
		"ex_rxn_C", "C", 1, "character",
		"pts_pd", "p", 2, "numeric",
		"output_root", "o", 2, "character",
		"skip_3d_compute", "s", 2, "logical",
		"skip_fva_compute", "f", 2, "logical",
		"ignore_trios_bounds", "i", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

PTS_PD=15;

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <model name>\n",
		"	-d <default exchanges, tsv>\n",
		"	-A <Exchange Reaction A>\n",
		"	-B <Exchange Reaction B>\n",
		"	-C <Exchange Reaction C>\n",
		"	[-p <points per direction, default=", PTS_PD, ">]\n",
		"	[-o <output root>]\n",
		"	[-s (skip 3D compute)]\n",
		"	[-f (skip FVA compute)]\n",
		"	[-i (ignore trios bounds and go from -1000 to 1000)]\n",
		"\n",
		"Script will read in the model and the 'default' exchange reaction bounds.\n",
		"Then for the three exchange reactions specified, combinations of the\n",
		"3 reactions fluxes at difference flux settings will be set and the objective\n",
		"function will be optimized.\n",
		"\n",
		"For example, if the lower and upper bound of the exchange reaction for A\n",
		"is set from -100 to 500, then the flux of reaction A will be set to combinations\n",
		"with exchange reaction for B and C, with A between -100 and 500.\n",
		"\n",
		"The points per direction (ppd) option will specify how many points should\n",
		"be computed in each direction, such that the total number of optimizations will\n",
		"be ppd^3.",
		"\n");

if(!length(opt$model) ||
	!length(opt$def_flux) ||
	!length(opt$ex_rxn_A) ||
	!length(opt$ex_rxn_B) ||
	!length(opt$ex_rxn_C) 
){
	cat(usage);
	q(status=-1);
}

library(sybil);
library(glpkAPI);
cat("\n\n");

Model=opt$model;
DefFlux=opt$def_flux;
ExRxnA=opt$ex_rxn_A;
ExRxnB=opt$ex_rxn_B;
ExRxnC=opt$ex_rxn_C;

Skip3DCompute=ifelse(length(opt$skip_3d_compute), T, F);
SkipFVACompute=ifelse(length(opt$skip_fva_compute), T, F);
IgnoreTriosBounds=ifelse(length(opt$ignore_trios_bounds), T, F);

PtsPD=ifelse(length(opt$pts_pd), opt$pts_pd, PTS_PD);
OutputRoot=ifelse(length(opt$output_root), opt$output_root, gsub("\\.tsv$", "", DefFlux));

cat("DefFlux File: ", DefFlux, "\n", sep="");
cat("Ex Rxn A: ", ExRxnA, "\n", sep="");
cat("Ex Rxn B: ", ExRxnB, "\n", sep="");
cat("Ex Rxn C: ", ExRxnC, "\n", sep="");
cat("\n");
cat("Points per Dimension: ", PtsPD, "\n", sep="");
cat("Output Root: ", OutputRoot, "\n", sep="");
cat("\n");

if(Skip3DCompute){
	cat("Skipping 3D Compute.\n");
}

# Clean the name just in case more than the prefix was specified.
Model=gsub("_[a-z]+\\.tsv$","",Model);
cat("Model: ", Model, "\n", sep="");

###############################################################################

read_exchange_rxns_bounds=function(fname){
	data=as.matrix(read.table(fname, sep="\t", header=T));
	bounds=as.data.frame(matrix(0, nrow=nrow(data), ncol=3));
	colnames(bounds)=c("lowbnd", "uppbnd", "name");
	rownames(bounds)=data[,1]
	bounds[,"lowbnd"]=as.numeric(data[,"lowbnd"]);
	bounds[,"uppbnd"]=as.numeric(data[,"uppbnd"]);
	bounds[,"name"]=data[,"name"];
	return(bounds);
}

apply_bounds_to_model=function(model, bounds){

	# bounds: 
	#	num_bnd_rxns x ("lowbnd", "uppbnd") matrix
	# 	row names are the rxn ids

	curS_lb=lowbnd(model);
	curS_ub=uppbnd(model);

	model_rxn_ids=react_id(model);
	bounds_rxn_ids=rownames(bounds);

	names(curS_lb)=model_rxn_ids;
	names(curS_ub)=model_rxn_ids;

	curS_lb[bounds_rxn_ids]=bounds[bounds_rxn_ids,"lowbnd"];
	curS_ub[bounds_rxn_ids]=bounds[bounds_rxn_ids,"uppbnd"];
	
	lowbnd(model)=curS_lb;
	uppbnd(model)=curS_ub;
	
	return(model);
}

test_model=function(model){
	opt_res=optimizeProb(model, algorithm="fba", retOptSol=TRUE);

	objective_flux=attributes(opt_res)$lp_obj;
	reaction_flux=as.vector(attributes(attributes(opt_res)$fluxdist)$fluxes);

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

cat("  Model Path: ", model_path, "\n");
cat("Model Prefix: ", model_prefix, "\n");
cat("\n");

#  model=readTSVmod(fpath="/local/devel/DAS/users/kli/SVN/DAS/SystemsBioTools/Plot_ObjectiveFlux_Ternary/TD", prefix="Thiobacillus_Denitrificans", suffix="tsv")

model=readTSVmod(fpath=model_path, prefix=model_prefix, suffix="tsv");
S=S(model);
Sdim=dim(S);

num_rxns=length(react_id(model));
num_mets=length(met_id(model));

cat("S Matrix dimensions: Row: ", Sdim[1], " Col: ", Sdim[2], "\n");
cat("Num Mets: ", num_mets, "\n");
cat("Num Rxns: ", num_rxns, "\n");
cat("\n");
cat("Example Reaction IDs:\n");
print(head(react_id(model)));
cat("Example Metabolite IDs:\n");
print(head(met_id(model)));
cat("\n");

###############################################################################

cat("Loading Default Exchange Reaction Bounds:\n");
def_flux_table=read_exchange_rxns_bounds(DefFlux);

# Generate a mapping from exchange reaction ID to exchange reaction description
id_to_desc_ex_rxn_map=as.matrix(def_flux_table)[,"name"];
names(id_to_desc_ex_rxn_map)=rownames(def_flux_table);

###############################################################################

tern_exrxns=c(ExRxnA, ExRxnB, ExRxnC);
tern_exrxn_bounds=def_flux_table[tern_exrxns, ];
cat("Exchange Reaction Bounds for Ternary Analysis:\n");
print(tern_exrxn_bounds);

cat("Applying default bounds to model:\n");
model_wdefbds=apply_bounds_to_model(model, def_flux_table);

# Allocate matrix for storing bounds
make_temp_flux_mat=function(){
	tmp_ex_flux_mat=matrix(0, nrow=3, ncol=2);
	colnames(tmp_ex_flux_mat)=c("lowbnd", "uppbnd");
	rownames(tmp_ex_flux_mat)=tern_exrxns;
	return(tmp_ex_flux_mat);
}

tmp_ex_flux_mat=make_temp_flux_mat();



###############################################################################

ZERO_TOL=1e-6;


find_zero_bounds=function(target_rxn_ID, optim_flux, model, search_bounds, max_trials){
	
	cat("\n");
	cat("Target Reaction ID: ", paste(target_rxn_ID, collapse=", "), "\n");
	cat("\tOptimum Trio Fluxes: ", paste(optim_flux, collapse=", "), "\n");
	trio_rxn_names=names(optim_flux);
	cat("\tTrio names: ", paste(trio_rxn_names, collapse=", "), "\n");

	for(direction in c(-1, 1)){

		cat("     Direction to search in: ", direction, "\n");

		if(direction==1){
			splits=seq(as.numeric(optim_flux[target_rxn_ID]), 
				as.numeric(search_bounds[2]), length.out=(max_trials+1));
		}else{
			splits=seq(as.numeric(optim_flux[target_rxn_ID]), 
				as.numeric(search_bounds[1]), length.out=(max_trials+1));
		}
		splits=splits[2:(max_trials+1)];
		
		cat("     Trial positions: ", paste(splits, collapse=", "), "\n");

		biomass_flux=1;
		num_trials=1;

		while(biomass_flux > ZERO_TOL && num_trials<=max_trials){

			# Set test value via bounds
			temp_flux_mat=make_temp_flux_mat();
			temp_flux_mat[trio_rxn_names, c("lowbnd","uppbnd")]=optim_flux[trio_rxn_names];
			temp_flux_mat[target_rxn_ID, c("lowbnd", "uppbnd")]=splits[num_trials];
			#print(temp_flux_mat);

			# Test model
			temp_model=apply_bounds_to_model(model, temp_flux_mat);
			test_flux=test_model(temp_model);
			biomass_flux=test_flux$objective_flux;
			cat("\tTest Flux Result: ", biomass_flux, "\n");

			num_trials=num_trials+1;
		}

		if(direction==1){
			ub_no_flux=splits[num_trials-1];
		}else{
			lb_no_flux=splits[num_trials-1];
		}
		
		
	}	

	cat("  Recommended Bounds: ", lb_no_flux, " to ", ub_no_flux, "\n");
	return(c(lb_no_flux, ub_no_flux));
	
}

###############################################################################

EstimateBounds=T;
QuickRefine=F;
ExpandBounds=F;

if(EstimateBounds){

	if(IgnoreTriosBounds){
		tmp_ex_flux_mat[ExRxnA,]=c(-1000,1000);
		tmp_ex_flux_mat[ExRxnB,]=c(-1000,1000);
		tmp_ex_flux_mat[ExRxnC,]=c(-1000,1000);
	}else{
		for(rxn_id in c(ExRxnA, ExRxnB, ExRxnC)){
			tmp_ex_flux_mat[rxn_id,]=c(tern_exrxn_bounds[rxn_id, "lowbnd"], tern_exrxn_bounds[rxn_id, "uppbnd"]);
		}
	}	

	unrestricted_model=apply_bounds_to_model(model_wdefbds, tmp_ex_flux_mat);
	optim_flux=test_model(unrestricted_model);
	cat("Optimum Biomass Flux: ", optim_flux$objective_flux, "\n");

	# Compute Flux at optimum
	optimum_flux_values=as.vector(optim_flux$reaction_flux);
	names(optimum_flux_values)=react_id(model);
	optim_flux_for_trio=optimum_flux_values[c(ExRxnA, ExRxnB, ExRxnC)];

	if(IgnoreTriosBounds){
		zero_bounds_A=c(-1000,1000);
		zero_bounds_B=c(-1000,1000);
		zero_bounds_C=c(-1000,1000);
	}else{
		zero_bounds_A=tern_exrxn_bounds[ExRxnA,];
		zero_bounds_B=tern_exrxn_bounds[ExRxnB,];
		zero_bounds_C=tern_exrxn_bounds[ExRxnC,];
	}

	max_trial_arr=c(4,5,5,10,10);

	if(QuickRefine){
		max_trial_arr=c(4,4);
	}

	for(i in 1:length(max_trial_arr)){
	# For each of member of the trio, try to refine bounds
		zero_bounds_A=find_zero_bounds(ExRxnA, optim_flux_for_trio, model_wdefbds, 
			zero_bounds_A, max_trial_arr[i]);
		zero_bounds_B=find_zero_bounds(ExRxnB, optim_flux_for_trio, model_wdefbds,
			zero_bounds_B, max_trial_arr[i]);
		zero_bounds_C=find_zero_bounds(ExRxnC, optim_flux_for_trio, model_wdefbds,
			zero_bounds_C, max_trial_arr[i]);
	}

	cat("\n");
	cat("Bounds for ", ExRxnA, " : ", zero_bounds_A[1], " to ", zero_bounds_A[2], "\n", sep="");
	cat("Bounds for ", ExRxnB, " : ", zero_bounds_B[1], " to ", zero_bounds_B[2], "\n", sep="");
	cat("Bounds for ", ExRxnC, " : ", zero_bounds_C[1], " to ", zero_bounds_C[2], "\n", sep="");

	expand_bounds=function(bounds, extra=0.2, maxflux=1000, minflux=-1000){
		lb=bounds[1];
		ub=bounds[2];
		cat("Before Extending: ", lb, " to ", ub , "\n");
		range=ub-lb;
		nom_extra=extra*range;
		lb=max(lb-nom_extra, minflux);
		ub=min(ub+nom_extra, maxflux);		
		cat("After Extending: ", lb, " to ", ub , "\n");
		return(c(lb,ub));	
	}

	if(ExpandBounds){
		zero_bounds_A=expand_bounds(zero_bounds_A);
		zero_bounds_B=expand_bounds(zero_bounds_B);
		zero_bounds_C=expand_bounds(zero_bounds_C);
	}

	incrA=seq(zero_bounds_A[1], zero_bounds_A[2], length.out=PtsPD);
	incrB=seq(zero_bounds_B[1], zero_bounds_B[2], length.out=PtsPD);
	incrC=seq(zero_bounds_C[1], zero_bounds_C[2], length.out=PtsPD);
}else{

	incrA=seq(tern_exrxn_bounds[ExRxnA, "lowbnd"], tern_exrxn_bounds[ExRxnA, "uppbnd"], length.out=PtsPD);
	incrB=seq(tern_exrxn_bounds[ExRxnB, "lowbnd"], tern_exrxn_bounds[ExRxnB, "uppbnd"], length.out=PtsPD);
	incrC=seq(tern_exrxn_bounds[ExRxnC, "lowbnd"], tern_exrxn_bounds[ExRxnC, "uppbnd"], length.out=PtsPD);
}

# Force inclusion of 0 point if the resultant bounds overlap 0.

add_zero=function(values){
	if(min(values)<0 && max(values)>0){
		cat("Added 0 to range.\n");
		values=sort(c(values,0));
	}
	return(values);
}

incrA=add_zero(incrA);
incrB=add_zero(incrB);
incrC=add_zero(incrC);

if(!Skip3DCompute){

	#Test that we can write the file before starting the computes
	write.table("Test to make sure file is writeable before starting computes.", 
		file=paste(OutputRoot, ".trios_flux.tsv", sep=""),
		row.names=F, col.names=F
		);

	tot=length(incrA)*length(incrB)*length(incrC);
	flux_val_mat=matrix(NA, nrow=tot, ncol=4);
	colnames(flux_val_mat)=c(ExRxnA, ExRxnB, ExRxnC, "OptFlux");

	ix=1;
	for(ai in incrA){
		for(bi in incrB){
			for(ci in incrC){
				tmp_ex_flux_mat[ExRxnA, ]=c(ai, ai);
				tmp_ex_flux_mat[ExRxnB, ]=c(bi, bi);
				tmp_ex_flux_mat[ExRxnC, ]=c(ci, ci);

				cat("\n[", ix, "/", tot, "] Working on:\n");
				print(tmp_ex_flux_mat);
				restricted_model=apply_bounds_to_model(model_wdefbds, tmp_ex_flux_mat);
				optim_flux=test_model(restricted_model);
				cat("Flux: ", optim_flux$objective_flux, "\n");

				flux_val_mat[ix, ]=c(ai, bi, ci, optim_flux$objective_flux);
				
				ix=ix+1;
			}
		}
	}

	write.table(flux_val_mat, file=paste(OutputRoot, ".trios_flux.tsv", sep=""), row.names=F, sep="\t", quote=F);

}

run_fva=function(model, target_rxn_ids){
# Runs the FVA and returns it in a matrix with two columns.
# target_rxn_ids, are indices
	fva_out=fluxVar(model, target_rxn_ids);
	num_indiv_exch_rxns=length(target_rxn_ids);
	fva_bounds=matrix(0, nrow=num_indiv_exch_rxns, ncol=2);
	colnames(fva_bounds)=c("lowbnd", "uppbnd");
        fva_bounds[,"lowbnd"]=attr(fva_out, "lp_obj")[1:num_indiv_exch_rxns];
        fva_bounds[,"uppbnd"]=attr(fva_out, "lp_obj")[(num_indiv_exch_rxns+1):(2*num_indiv_exch_rxns)];
	return(fva_bounds);
}

get_indices_for_rxns_by_rxn_id=function(model, targeted_rxn_ids){
	rxn_ids=react_id(model);
	num_rxn_ids=length(targeted_rxn_ids);

	indices=numeric(num_rxn_ids);	
	errors=0;
	for(i in 1:num_rxn_ids){
		matching=targeted_rxn_ids[i]==rxn_ids;
		idx=which(matching);
		if(length(idx)==0){
			cat("Could not find: ", targeted_rxn_ids[i], " in model.\n");
			errors=errors+1;
		}else{
			indices[i]=which(targeted_rxn_ids[i]==rxn_ids);	
		}
	}

	if(errors>0){
		cat("Error: could not find matching IDs for reactions.\n");
		quit(status=-1);
	}

	return(indices);
}


plot_array_of_bounds=function(lb_arr, ub_arr, col=c("black"), lblim, ublim, 
	label_flux=F, label_xaxes=F, idmap=NULL, target_trio=""){

	num_pairs=length(lb_arr);
	left_right_margin=.75;

	if(num_pairs>10){
		num_guide_lines=7;
		guide_lines=seq(.5, num_pairs+.5, 1);
		guide_idx=round(seq(1, length(guide_lines), length.out=num_guide_lines+2)); 	
		guide_lines=guide_lines[guide_idx];
		guide_lines=guide_lines[-c(1,num_guide_lines+2)];
	}else{
		guide_lines=c();
	}

	rxn_name=names(lb_arr);	

	if(!label_xaxes){
		plot(0, type="n", ylim=c(lblim, ublim), xlim=c(1-left_right_margin, num_pairs+left_right_margin),
			 yaxt="n", xaxt="n");
		abline(h=0, col="grey");
		width=.75;
		half_width=width/2;
		if(length(col)==1){
			col=rep(col, num_pairs);
		}
		for(i in 1:num_pairs){

			# If rxn_name was targeted for analysis, highlight it's position
			if(rxn_name[i]==target_trio){
				abline(v=i, lwd=16, col="khaki");
			}

			rect(xleft=i-half_width, xright=i+half_width, 
				ybottom=0, ytop=ub_arr[i], col="lightblue3",
				border="lightblue3");
			rect(xleft=i-half_width, xright=i+half_width, 
				ybottom=lb_arr[i], ytop=ub_arr[i], col=col[i],
				border=col[i]);


			if(label_flux){
				if(lb_arr[i]>=0){
					text(i, 0, label=sprintf("%4.4g", lb_arr[i]), srt=-90, adj=c(-.1,.25), cex=.75);
				}else{
					text(i, 0, label=sprintf("%4.4g", lb_arr[i]), srt=90, adj=c(-.1,.35), cex=.75);
				}
			}
		}
	}else{
		# Plot the names of each reaction
		plot(0, type="n", ylim=c(0, 1000), xlim=c(1-left_right_margin, num_pairs+left_right_margin),
			 yaxt="n", xaxt="n");
		size=min(1, 27/num_pairs);
		for(i in 1:num_pairs){
			label_name=idmap[rxn_name[i]]
			if(is.na(label_name)){
				label_name=rxn_name[i];
			}
			text(x=i, y=1000, label=label_name, srt=-90, pos=4, offset=0.0, cex=size);
		}
	}
	abline(v=guide_lines, col="grey", lty="dotted");
}

plot_fva=function(trios_settings, flux, posneg_exchanges, pos_only_exchanges, label_xaxes=F, idmap=NULL, target_trio=""){
	#print(trios_settings);
	#print(flux);
	#print(posneg_exchanges);
	#print(pos_only_exchanges);

	# Plot trios and biomass
	cat("Plotting Trios and Objective...\n");
	trio_and_obj=c(trios_settings, flux);
	names(trio_and_obj)=c(names(trios_settings), "Objective");
	plot_array_of_bounds(trio_and_obj, c(0,0,0,0,0), lblim=-1100, ublim=1100, 
		col=c("red","green","blue", "black"), label_flux=T, label_xaxes=label_xaxes, idmap=idmap, 
		target_trio=target_trio);


	num_posneg=nrow(posneg_exchanges);
	num_pos_only=nrow(pos_only_exchanges);

	if(flux<ZERO_TOL && label_xaxes==F){
		posneg_exchanges[,"lowbnd"]=rep(0, num_posneg);
		posneg_exchanges[,"uppbnd"]=rep(0, num_posneg);
		pos_only_exchanges[,"lowbnd"]=rep(0, num_pos_only);
		pos_only_exchanges[,"uppbnd"]=rep(0, num_pos_only);
	}

	# Plot pos/neg exchanges
	cat("Plotting Pos/Neg Exchanges...\n");
	plot_array_of_bounds(posneg_exchanges[,"lowbnd"], posneg_exchanges[,"uppbnd"], lblim=-1100, ublim=1100,
		label_flux=F, label_xaxes=label_xaxes, idmap=idmap);


	if(!label_xaxes){
		marks=c(-1000,-500, 500, 1000);
		axis(side=2, at=marks, labels=marks, cex.axis=.4, line=-2, tick=F, las=2);
	}

	# Plot positive only bounds
	cat("Plotting Pos Only Exchanges...\n");
	plot_array_of_bounds(pos_only_exchanges[,"lowbnd"], pos_only_exchanges[,"uppbnd"], lblim=-100, ublim=1100,
		label_flux=F, label_xaxes=label_xaxes, idmap=idmap);

	if(!label_xaxes){
		marks=c(250, 500, 750, 1000);
		axis(side=2, at=marks, labels=marks, cex.axis=.4, line=-2, tick=F, las=2);
	}

}

build_layout=function(num_rows){
	# Constructs a layout based on the number of rows assuming the 3 horizontal graphs
	layout_line=c(1,2,2,2,3,3,3);
	line_len=length(layout_line);
	layout_matrix=matrix(0, nrow=num_rows, ncol=line_len);
	for(i in 1:num_rows){
		layout_matrix[i,]=layout_line+(i-1)*3;
	}
	#print(layout_matrix);
	return(layout_matrix);	
}

match_order=function(src_matrix, preferred_order_ids){

	src_names=rownames(src_matrix);
	src_names_len=length(src_names);
	num_preferred_set=length(preferred_order_ids);
	
	kept_names_in_order=character(src_names_len);
	
	shared_ix=1;
	for(i in 1:num_preferred_set){
		cur_id=preferred_order_ids[i];
		if(length(intersect(cur_id, src_names))){
			kept_names_in_order[shared_ix]=cur_id;	
			shared_ix=shared_ix+1;
		}
	}

	return(src_matrix[kept_names_in_order,]);
}



if(!SkipFVACompute){
	# FVA compute at optimal axes
	cat("\n");
	cat("Optimal Flux for Trio:\n"); print(optim_flux_for_trio);
	cat("\n");
	cat("Ranges for FVA analysis: \n"); cat(ExRxnA, ":\n");
	print(incrA);
	cat(ExRxnB, ":\n");
	print(incrB);
	cat(ExRxnC, ":\n");
	print(incrC);
	cat("\n");

	trio_rxn_ids=c(ExRxnA, ExRxnB, ExRxnC);

	# Convert default flux data.frame into matrix so we can access it by named idx
	def_flux_matrix=as.matrix(def_flux_table[,c("lowbnd","uppbnd")]);

	# Remove trios IDs from set we will perform FVA on
	all_def_flux_names=rownames(def_flux_matrix);
	def_flux_names_less_trio=setdiff(all_def_flux_names, trio_rxn_ids);
	def_flux_matrix=def_flux_matrix[def_flux_names_less_trio,];
	def_flux_table_names=rownames(def_flux_matrix);

	# Split exchange fluxes by what we ranges
	pos_flux_only=(def_flux_matrix[,"lowbnd"]>=0) & (def_flux_matrix[,"uppbnd"]>=0);
	posneg_flux=(def_flux_matrix[,"lowbnd"]<0) & def_flux_matrix[,"lowbnd"]<def_flux_matrix[,"uppbnd"];

	pos_flux_only_rxn_names=def_flux_table_names[pos_flux_only];
	posneg_flux_rxn_names=def_flux_table_names[posneg_flux];

	cat("Exchanges with only positive flux: \n");
	print(pos_flux_only_rxn_names);
	pos_flux_only_rxn_indices=get_indices_for_rxns_by_rxn_id(model_wdefbds, pos_flux_only_rxn_names);
	#print(pos_flux_only_rxn_indices);
	cat("\n");
	cat("Exchanges with positive or negative flux: \n");
	print(posneg_flux_rxn_names);
	posneg_flux_rxn_indices=get_indices_for_rxns_by_rxn_id(model_wdefbds, posneg_flux_rxn_names);
	#print(posneg_flux_rxn_iddndices);
	cat("\n");

	rxn_flux_pts=list();
	rxn_flux_pts[[ExRxnA]]=incrA;
	rxn_flux_pts[[ExRxnB]]=incrB;
	rxn_flux_pts[[ExRxnC]]=incrC;

	pdf(paste(OutputRoot, ".trios_fva.pdf", sep=""), height=11, width=8.5);
	par(oma=c(0,0,2,0));
	
	for(cur_rxn_id in trio_rxn_ids){
		cat("\n\n");
		cat("Performing FVA on ",  cur_rxn_id, "\n", sep="");
		cur_rxn_pts=rxn_flux_pts[[cur_rxn_id]];

		num_rxn_pts=length(cur_rxn_pts);
		par(mar=c(0,0,0,0));
		layout(build_layout(num_rxn_pts+1));

		for(i in 1:num_rxn_pts){
			
			# Make and populate trios settings in data structure
			cur_ex_flux_mat=make_temp_flux_mat();
			cur_ex_flux_mat[trio_rxn_ids, c("lowbnd","uppbnd")]=optim_flux_for_trio[trio_rxn_ids];
			cur_ex_flux_mat[cur_rxn_id, c("lowbnd", "uppbnd") ]=cur_rxn_pts[i];

			# Apply trios settings to model
			print(cur_ex_flux_mat);
			restricted_model=apply_bounds_to_model(model_wdefbds, cur_ex_flux_mat);

			# Compute optimal flux
			optim_flux=test_model(restricted_model);
			cat("Flux: ", optim_flux$objective_flux, "\n");

			# Compute FVA for pos flux and pos/neg flux separately
			pos_flux_only_fva=run_fva(restricted_model, pos_flux_only_rxn_indices );
			rownames(pos_flux_only_fva)=pos_flux_only_rxn_names;
			posneg_flux_fva=run_fva(restricted_model, posneg_flux_rxn_indices );
			rownames(posneg_flux_fva)=posneg_flux_rxn_names;

			# Assemble trios values used in iteration
			trio_setting=optim_flux_for_trio[trio_rxn_ids];
			trio_setting[cur_rxn_id]=cur_rxn_pts[i];

			# Sort flux orders by input
			posneg_flux_fva=match_order(posneg_flux_fva, def_flux_table_names);
			pos_flux_only_fva=match_order(pos_flux_only_fva, def_flux_table_names);

			# Plot FVA for this iteration
			plot_fva(trio_setting, optim_flux$objective_flux, 
				posneg_flux_fva, pos_flux_only_fva, label_xaxes=F, target_trio=cur_rxn_id);
		}	

		# Plot the labels on the bottom of page
		plot_fva(trio_setting, optim_flux$objective_flux,
			 posneg_flux_fva, pos_flux_only_fva, label_xaxes=T, idmap=id_to_desc_ex_rxn_map);

		# Label page in top margin
		mtext(paste(OutputRoot, " (", cur_rxn_id, ") ", sep=""), outer=T);
	}

}

###############################################################################

print(warnings());
cat("Done.\n");


























