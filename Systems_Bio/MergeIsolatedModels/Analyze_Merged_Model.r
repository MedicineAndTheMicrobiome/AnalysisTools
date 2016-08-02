#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
		"model", "m", 1, "character",
		"biomass_rxn_A", "A", 1, "character",
		"biomass_rxn_B", "B", 1, "character",
		"name_A", "a", 1, "character",
		"name_B", "b", 1, "character",
		"output", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <model name>\n",
		"	-A <biomass reaction ID of model A>\n",
		"	-B <biomass reaction ID of model B>\n",
		"	-a <name of model A>\n",
		"	-b <name of model B>\n",
		"	[-o <output filename root>]\n",
		"\n",
		"This script will read in single merged model then\n",
		"perform optimizations across a range of biomass ratios\n",
		"by changing the proportions of stoichiometries in\n",
		"each of the individual model's biomass reaction.\n",
		"\n",
		"The -m specifies the model which is a TSV file with reactions\n",
		"equations, flux bounds, etc., in it.\n",
		"\n",
		"The -A and -B parameter should specify which reaction ID should\n",
		"be considered the biomass reaction for each organism\n",
		"respectively.\n",
		"\n",
		"The -a and -b parameter is used for labeling in the\n",
		"generated graphs.\n",
		"\n",
		"The following output files will be generated:\n",
		"	1.)  comb.rxn_flux.tsv\n",
		"		Contains the fluxes across all proportions for a\n",
		"		single possible solution.\n",
		"	2.)  comb.fva.pdf\n",
		"		Contains plot of combined biomass across all proportions\n",
		"		then plots of FVA across shared exchanges\n",
		"	3.)  comb.opt_flux.tsv\n",
		"		Contains the values for the optimal flux across\n",
		"		across all proportions.\n",
		"	4.)  indiv.fva.pdf\n",
		"		Similar to the comb.fva.pdf, except with the addition that\n",
		"		the flux variability across each individual organism is plotted\n",
		"		on each side of the combined flux in different colors.\n",
		"\n");

if(!length(opt$model) ){
	cat(usage);
	q(status=-1);
}

library(sybil);
library(glpkAPI);
cat("\n\n");

Model_name=opt$model;
BiomassID_A=opt$biomass_rxn_A;
BiomassID_B=opt$biomass_rxn_B;
NameA=opt$name_A;
NameB=opt$name_B;

Model_name=gsub("_[a-z]+\\.tsv$","",Model_name);

if(length(opt$output)>0){
	Output=opt$output;
}else{
	Output=Model_name;
}

# Clean the name just in case more than the prefix was specified.
cat("Model: ", Model_name, "\n", sep="");
cat("\n");
cat("Biomass A: ", BiomassID_A, "\n", sep="");
cat("Name A: ", NameA, "\n", sep="");
cat("\n");
cat("Biomass B: ", BiomassID_B, "\n", sep="");
cat("Name B: ", NameB, "\n", sep="");
cat("\n");
cat("Output Filename Root: ", Output, sep="");
cat("\n");

###############################################################################
# Load model

# Find the model
path_arr=strsplit(Model_name, "/")[[1]];
model_prefix=tail(path_arr, 1);
model_path=paste(head(path_arr, length(path_arr)-1), collapse="/");
if(model_path==""){
	model_path=paste(".",model_path, sep="");
}

# Get file info on model
met_info=file.info(paste(Model_name, "_met.tsv", sep=""));
react_info=file.info(paste(Model_name, "_react.tsv", sep=""));
desc_info=file.info(paste(Model_name, "_desc.tsv", sep=""));

cat("\n");
cat("Last Modified:\n");
cat("  Description: ", format(as.POSIXct(desc_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"), "\n");
cat("  Metabolites: ", format(as.POSIXct(met_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"),  "\n");
cat("    Reactions: ", format(as.POSIXct(react_info$mtime, origin="1970-01-01"), format="%B %d %Y, %H:%M:%S"), "\n");
cat("\n");

cat("  Model Path: ", model_path, "\n");
cat("Model Prefix: ", model_prefix, "\n");
cat("\n");

###############################################################################

test_model=function(model){
	opt_res=optimizeProb(model, algorithm="fba", retOptSol=TRUE);
	#print(opt_res);

	objective_flux=attributes(opt_res)$lp_obj;
	reaction_flux=(attributes(attributes(opt_res)$fluxdist)$fluxes);

	optim_out=list();
	optim_out$objective_flux=objective_flux;
	optim_out$reaction_flux=reaction_flux;
	
	return(optim_out);
}

###############################################################################

# Load the model
model=readTSVmod(fpath=model_path, prefix=model_prefix, suffix="tsv");
original_model=model;

# model=readTSVmod(fpath=".", prefix="Combined", suffix="tsv"); original_model=model;
# BiomassID_A="BIO_Rfer3.RFe"; BiomassID_B="ecoli_biomass.TDe";
# NameA="Rhodoferax_Ferrireducens"; NameB="Thiobacillus_Denitrificans";


S=S(model); #S[metabolites, reactions];
react_id=react_id(model);
obj_coefs=obj_coef(model);

biomassA_idx=which(react_id==BiomassID_A);
biomassB_idx=which(react_id==BiomassID_B);

A_stoich=S[,biomassA_idx];
B_stoich=S[,biomassB_idx];

# Turn off A as objective function and reaction
obj_coefs[biomassA_idx]=0; obj_coef(model)=obj_coefs;
lowbnds=lowbnd(model); lowbnds[biomassA_idx]=0; lowbnd(model)=lowbnds;
uppbnds=uppbnd(model); uppbnds[biomassA_idx]=0; uppbnd(model)=uppbnds;

# Store combined stoichiometry in B
increments=seq(0,1,.05);
#increments=seq(0,1,.5);
num_increments=length(increments);

combined_optimal_flux=numeric(num_increments);

find_exchange_reactions=function(Smat){
	ins=(Smat<0);
	outs=(Smat>0);

	num_reactions=ncol(Smat);
	any_ins=apply(ins, 2, sum)>0;
	any_outs=apply(outs, 2, sum)>0;

	exch_rxns_idx=which(any_ins & (any_outs==F));
	return(exch_rxns_idx);
}

exch_rxn_idx=find_exchange_reactions(S);
cat("Exchange Reaction Indices: \n");
print(exch_rxn_idx);
num_fva_rxns=length(exch_rxn_idx);

ZERO_THRESHOLD=0.0001;

lb_mat=array(dim=c(num_fva_rxns, num_increments));
ub_mat=array(dim=c(num_fva_rxns, num_increments));

num_reactions=ncol(S(model));
rxn_flux_mat=matrix(0, nrow=num_reactions, ncol=num_increments);

for(i in 1:num_increments){

	# Prepare new stoichiometry matrix with weighted biomass proportions
	p=increments[i];
	cat("Working on ", sprintf("%03.02g", p) , " ", NameA, " vs. ", sprintf("%03.02g",(1-p)), " ", NameB, "\n", sep="");
	combined_stoich=A_stoich*p + B_stoich*(1-p);
	S[,biomassB_idx]=combined_stoich;
	S(model)=S;
	
	# Run and store FBA
	optim_out=test_model(model);
	combined_optimal_flux[i]=optim_out$objective_flux;
	rxn_flux_mat[,i]=optim_out$reaction_flux[,1];

	cat("Biomass Flux: ", rxn_flux_mat[biomassB_idx, i], "\n");

	# Run and store FVA
	fva_out=fluxVar(model, exch_rxn_idx);
	lb_mat[,i]=attr(fva_out, "lp_obj")[1:num_fva_rxns];
	ub_mat[,i]=attr(fva_out, "lp_obj")[(num_fva_rxns+1):(2*num_fva_rxns)];
}

generate_reaction_string=function(stoich, met_ids, flux_dir){
	
	lhs=which(stoich <= -1);
	rhs=which(stoich >= 1);

	lhs_coef=stoich[lhs]*-1;
	rhs_coef=stoich[rhs];

	lhs_str=paste("(", lhs_coef, ")",  met_ids[lhs], collapse=" + ", sep="");
	rhs_str=paste("(", rhs_coef, ")",  met_ids[rhs], collapse=" + ", sep="");

	if(flux_dir==1){
		rxn_str=paste(lhs_str, " --> ", rhs_str, sep="");
	}else if(flux_dir==-1){
		rxn_str=paste(rhs_str, " --> ", lhs_str, sep="");
	}else{
		rxn_str=paste(lhs_str, " <=> ", rhs_str, sep="");
	}

	rxn_str=gsub("\\(1\\)","", rxn_str);

	return(rxn_str);
}

output_reactions=function(flux_mat, model, fname){

	num_trials=ncol(flux_mat);
	num_rxns=nrow(flux_mat);

	# Get names and ids
	react_names=react_name(model);
	react_ids=react_id(model);
	met_ids=met_id(model);

	# Determine which reactions have any flux across all combinations
	flux_mat[abs(flux_mat)<ZERO_THRESHOLD]=0;
	any_nonzero_flux=apply(flux_mat!=0, 1, any)

	fh=file(fname, "w");
	
	# Print header
	cat(file=fh, paste(c(
		"Reaction_ID",
		"Reaction_Name",
		"Reaction",
		"Flipping"
		), collapse="\t"));

	for(i in 1:num_trials){
		cat(file=fh, "\t", sprintf("%03.02g%%", increments[i]*100));
	}
	cat(file=fh, "\n");

	# Print reaction fluxes across all trials
	for(i in 1:num_rxns){
		if(any_nonzero_flux[i]){

			# Determine if there is a consensus flux direction
			all_for=all(flux_mat[i,] >= 0);
			all_rev=all(flux_mat[i,] <= 0);
			flipping="--";
			if(all_for){
				flux_dir=1;
			}else if(all_rev){
				flux_dir=-1;
			}else{
				flux_dir=0;
				flipping="Flipping";
			}

			rxn_str=generate_reaction_string(S[,i], met_ids, flux_dir);

			cat(file=fh, paste(c(react_ids[i], react_names[i], rxn_str, flipping), collapse="\t"));
			cat(file=fh, "\t");
			cat(file=fh, paste(flux_mat[i,], collapse="\t"));
			cat(file=fh, "\n");
		}
	}

	close(fh);
}

flux_fname=paste(Output, ".comb.rxn_flux.tsv", sep="");
output_reactions(rxn_flux_mat, model, flux_fname);

###############################################################################

plot_flux=function(flux_arr, nameA, nameB, title=""){
	num_flux_values=length(flux_arr);
	cat("Num flux values: ", num_flux_values, "\n", sep="");

	max_flux=max(flux_arr);
	max_pos=which(max_flux==flux_arr);

	plot(flux_arr, main="Combined Biomass", type="b", xaxt="n", 
		xlim=c(0, num_flux_values+1),
		xlab="", ylab="Flux", ylim=c(0, max_flux*1.2));
	mtext(title, side=3, outer=T, cex=.8, adj=0);
	text(max_pos, max_flux, adj=c(-.15,-1), labels=sprintf("Max Flux: %5.3f", max_flux), col="blue");
	abline(v=max_pos, col="blue", lty="dashed");

	prop=(max_pos-1)/(num_flux_values-1);
	text(max_pos, 0, labels=sprintf("Max Ratio: %2.2g%% / %2.2g%%", prop*100, (1-prop)*100), col="blue", cex=.7, pos=4);
	
	axis(1, at=c(1, (num_flux_values+1)/2, num_flux_values), 
		labels=c(paste("100% ", nameB, sep=""),
			 "50%/50%", 
			paste("100% ", nameA, sep="")), cex=.4);

}

plot_fva=function(flux_lb, flux_ub, rxn_names, nameA, nameB){

	num_rxns=nrow(flux_lb);
	num_models=ncol(flux_lb);

	cat("Num Models: ", num_models, "\n", sep="");
	cat("Num Reactions: ", num_rxns, "\n", sep="");

	page_number=1;

	for(rxn_idx in 1:num_rxns){

		all_points=c(flux_lb[rxn_idx,], flux_ub[rxn_idx,])
		abs_max=max(abs(all_points));

		# Only plot exhanges that are non-zero
		if(any(abs(all_points)>ZERO_THRESHOLD)){

			# Page number
			mfg=par()$mfg;
			if(mfg[1]==mfg[3] && mfg[2]==mfg[4]){
				mtext(page_number, side=1, outer=T, cex=.8);
				page_number=page_number+1;
			}

			# Open a new plot with the ranges we want
			plot(0, 
				main=rxn_names[rxn_idx], type="n", xaxt="n", xlab="", ylab="Flux",
				xlim=c(0, num_models+1),
				ylim=c(-abs_max*1.1, abs_max*1.1)
			);

			# Lable the 0 flux line
			abline(h=0, col="grey", lty="dashed");

			# Label the x-axis
			axis(1, at=c(1, (num_models+1)/2, num_models), 
				labels=c(paste("100% ", nameB, sep=""),
					 "50%/50%", 
					paste("100% ", nameA, sep="")), cex=.4);

			# Plot a bar for each model
			for(model_ix in 1:num_models){
				
				lb=flux_lb[rxn_idx,model_ix]; ub=flux_ub[rxn_idx,model_ix];

				color="black";
				if(lb<=0 && ub<=0){ color="red"; }
				if(lb>=0 && ub>=0){ color="green"; }
				if(lb==0 && ub==0){ color="black"; }
				
				lines(
					c(model_ix, model_ix), 
					c(lb, ub),
					type="b", pch=0,
					col=color
				);
			}

		}
	}
	mtext(page_number, side=1, outer=T, cex=.8);
	
}

#------------------------------------------------------------------------------

pdf(paste(Output, ".comb.fva.pdf", sep=""),  height=11, width=8.5);
par(oma=c(1,0,1,0));
par(mfrow=c(4,1))
plot_flux(combined_optimal_flux, NameA, NameB, Output);
rxn_names=react_name(model)[exch_rxn_idx];
plot_fva(lb_mat, ub_mat, rxn_names, NameA, NameB);

dev.off();

#------------------------------------------------------------------------------

fh=file(paste(Output, ".comb.opt_flux.tsv", sep=""), "w");

max_flux=max(combined_optimal_flux);
max_pos=which(max_flux==combined_optimal_flux);
opt_prop=(max_pos-1)/(num_increments-1);
proportions=sprintf("%.2g%%", 100*(0:(num_increments-1))/(num_increments-1));

cat(file=fh, paste(c("# ModelName", "NameA", "NameB", "MaxFlux", "OptProp", proportions), collapse="\t"), "\n", sep="");
cat(file=fh, paste(c(Output, NameA, NameB, max_flux, opt_prop, combined_optimal_flux), collapse="\t"), "\n", sep="");

close(fh);

###############################################################################
###############################################################################

find_exchange_into_shared=function(shared_exch_rxn_idx, model, Aextcomp, Bextcomp){
	Smat=S(model);
	met_ids=met_id(model);
	rxn_names=react_name(model);
	met_comps=met_comp(model);
	mod_comparts=mod_compart(model);

	# Get compart id of A and B
	mod_compart_ix=c(which(mod_comparts==Aextcomp), which(mod_comparts==Bextcomp));
	cat("Compartment of A and B: ", paste(mod_compart_ix, collapse=" / "), "\n");

	# Get metabolites of interest
	share_exch_rxns=Smat[,shared_exch_rxn_idx];
	num_exch_rxns=ncol(share_exch_rxns);
	cat("Num exchange reactions: ", num_exch_rxns, "\n");

	# Fill in this matrix with the individual exchanges with the shared
	rxn_matrix=matrix(0, nrow=num_exch_rxns, ncol=4);
	colnames(rxn_matrix)=c("MetID", "A_Exchange", "B_Exchange", "Combined_Exchange");

	# Identify reactions that are only are channeling metabolites in/out
	passive_transport_rxns=apply(Smat, 2, sum)==0 & apply(abs(Smat), 2, sum)==2;

	for(i in 1:num_exch_rxns){
		dst_met_idx_of_shared_rxn=which(share_exch_rxns[,i]==-1);
		cat("[", i, "] Looking for: ", 
			met_ids[dst_met_idx_of_shared_rxn], 
			" (", dst_met_idx_of_shared_rxn, ") \n", sep="");

		# Get reactions involving shared reaction
		rxns_ix_with_shared_met=which(Smat[dst_met_idx_of_shared_rxn,]==1 & passive_transport_rxns);
		#cat("Rxn ID: ", rxns_ix_with_shared_met, "\n");
		#cat("Rxn ID: ", paste(rxn_names[rxns_ix_with_shared_met], collapse=" / "), "\n");

		# Get metabolite id of 
		src_met_idx=c();
		for(idx in rxns_ix_with_shared_met){
			src_met_idx=c(src_met_idx, which(Smat[,idx]==-1));
		}
		cat("Met ID: ", paste(met_ids[src_met_idx], collapse=" / "), "\n");

		# Get compartment of metabolite
		src_met_comp_idx=met_comps[src_met_idx];
		#cat("Comp ID: ", paste(src_met_comp_idx, collapse=" / "), "\n");

		rxn_matrix[i,1]=dst_met_idx_of_shared_rxn;

		if(!is.na(src_met_comp_idx[1])){
			if(mod_compart_ix[1]==src_met_comp_idx[1]){
				rxn_matrix[i, 2]=rxns_ix_with_shared_met[1];
			}else if(mod_compart_ix[2]==src_met_comp_idx[1]){
				rxn_matrix[i, 3]=rxns_ix_with_shared_met[1];
			}
		}

		if(!is.na(src_met_comp_idx[2])){
			if(mod_compart_ix[1]==src_met_comp_idx[2]){
				rxn_matrix[i, 2]=rxns_ix_with_shared_met[2];
			}else if(mod_compart_ix[2]==src_met_comp_idx[2]){
				rxn_matrix[i, 3]=rxns_ix_with_shared_met[2];
			}
		}

		rxn_matrix[i, 4]=shared_exch_rxn_idx[i];
		cat("\n");
	}
	return(rxn_matrix);
}

#------------------------------------------------------------------------------


plot_individual_exch_rxns=function(rxns_assoc_mat, rxn_ids_for_fva, flux_lb, flux_ub, model, nameA, nameB){

	met_names=met_name(model);
	num_shared_rxns=nrow(rxns_assoc_mat);
	num_models=ncol(flux_lb);
	lowbnds=lowbnd(model);
	uppbnds=uppbnd(model);

	cat("Num Models: ", num_models, "\n", sep="");
	cat("Num Shared Exchange Rxns: ", num_shared_rxns, "\n", sep="");

	#-----------------------------------------------------------------------------
	# Reorder
	preferred=c(
		"O2",
		"H+",
		"H2O",
		"Fe2+",
		"Fe3+",
		"pyrite FeS2",
		"Ammonium",
		"Nitrate",
		"Nitrite",
		"NO-c0",
		"Nitrous-oxide-c0",
		"N2-c0",
		"Urea",
		"Thiosulfate",
		"Sulfate",
		"Sulfite",
		"Hydrogen sulfide",
		"CO2",
		"Pyruvate",
		"Acetate",
		"Succinate",
		"Fumarate"
	);

	len_preferred=length(preferred);
	preferred_order=numeric(len_preferred);
	met_names_used=met_names[rxns_assoc_mat[,1]]

	for(i in 1:len_preferred){
		match_ix=which(met_names_used==preferred[i]);
		if(length(match_ix!=0)){
			preferred_order=c(preferred_order, match_ix);
		}
	}
	
	not_specified=(1:num_shared_rxns)[-preferred_order];
	preferred_order=c(preferred_order, not_specified);

	#-----------------------------------------------------------------------------

	offset=0.15;
	page_number=1;

	#for(rxn_idx in 1:num_shared_rxns){
	for(rxn_idx in preferred_order){

		shared_met_idx=rxns_assoc_mat[rxn_idx,1];
		flux_rxn_idx_A=rxns_assoc_mat[rxn_idx,2];
		flux_rxn_idx_B=rxns_assoc_mat[rxn_idx,3];
		flux_rxn_idx_comb=rxns_assoc_mat[rxn_idx,4];
		cat("Flux Rxn IDs: ", paste(rxns_assoc_mat[rxn_idx,], collapse=", "), "\n");

		comb_fva_idx=which(rxn_ids_for_fva==flux_rxn_idx_comb);
		A_fva_idx=which(rxn_ids_for_fva==flux_rxn_idx_A);
		B_fva_idx=which(rxn_ids_for_fva==flux_rxn_idx_B);

		comb_lb=flux_lb[comb_fva_idx,];
		A_lb=flux_lb[A_fva_idx,];
		B_lb=flux_lb[B_fva_idx,];

		comb_ub=flux_ub[comb_fva_idx,];
		A_ub=flux_ub[A_fva_idx,];
		B_ub=flux_ub[B_fva_idx,];

		model_constr_lb=lowbnds[flux_rxn_idx_comb];
		model_constr_ub=uppbnds[flux_rxn_idx_comb];

		all_points=c(comb_lb, A_lb, B_lb, comb_ub, A_ub, B_ub);
		abs_max=max(abs(all_points));

		# Only plot exhanges that are non-zero
		if(any(abs(all_points)>ZERO_THRESHOLD)){

			# Open a new plot with the ranges we want
			plot(0, main=met_names[shared_met_idx], type="n", xaxt="n", xlab="", ylab="Flux",
				xlim=c(0, num_models+1), ylim=c(-abs_max*1.1, abs_max*1.1));

			# Page number
			mfg=par()$mfg;
			print(mfg);
			if(mfg[1]==mfg[3] && mfg[2]==mfg[4]){
				mtext(page_number, side=1, outer=T, cex=.8);
				page_number=page_number+1;
			}

			# Lable the 0 flux line
			abline(h=0, col="grey", lty="dashed");

			# Label the x-axis
			#axis(1, at=c(1, (num_models+1)/2, num_models), 
			#	labels=c(paste("100% ", nameA, sep=""), "50%/50%", paste("100% ", nameB, sep="")),
			#	 col.axis=c("blue"), cex=.4);
			axis(1, at=c(1), 
				labels=paste("100% ", nameB, sep=""),
				 col.axis=c("blue"), cex=.4);
			axis(1, at=c((num_models+1)/2), 
				labels="50%/50%",
				 col.axis=c("black"), cex=.4);
			axis(1, at=c(num_models), 
				labels=paste("100% ", nameA, sep=""),
				 col.axis=c("dark green"), cex=.4);

			# Draw constraint bounds
			constr_col="red";
			# Lower bounds
			lines(c(0, num_models+1), c(model_constr_lb, model_constr_lb), col=constr_col, lty="dashed");
			points(c(0, num_models+1), c(model_constr_lb, model_constr_lb), 
				col=constr_col, bg=constr_col, pch=24, cex=1.2);
			# Upper bounds
			lines(c(0, num_models+1), c(model_constr_ub, model_constr_ub), col=constr_col, lty="dashed");
			points(c(0, num_models+1), c(model_constr_ub, model_constr_ub), 
				col=constr_col, bg=constr_col, pch=25, cex=1.2);

			# Plot a bar for each model
			for(model_ix in 1:num_models){

				# Plot combined flux
				lb=comb_lb[model_ix]; ub=comb_ub[model_ix];
				lines(c(model_ix, model_ix), c(lb, ub), 
					type="l", lwd=3, lend="butt", col="gray50");
				same=(abs(lb-ub) < abs_max*.01);
				points(model_ix[same], lb[same], pch="+", col="gray50", font=2);
				points(model_ix[same], ub[same], pch="+", col="gray50", font=2);

				# Plot A's flux
				lb=A_lb[model_ix]; ub=A_ub[model_ix];
				lines(c(model_ix+offset, model_ix+offset), c(lb, ub), 
					type="l", lwd=2, lend="butt", col="dark green");
				same=(abs(lb-ub) < abs_max*.01);
				points(model_ix[same]+offset, lb[same], pch="+", col="dark green", font=2);
				points(model_ix[same]+offset, ub[same], pch="+", col="dark green", font=2);

				# Plot B's flux
				lb=B_lb[model_ix]; ub=B_ub[model_ix];
				lines(c(model_ix-offset, model_ix-offset), c(lb, ub), 
					type="l", lwd=2, lend="butt", col="blue");
				same=(abs(lb-ub) < abs_max*.01);
				points(model_ix[same]-offset, lb[same], pch="+", col="blue", font=2);
				points(model_ix[same]-offset, ub[same], pch="+", col="blue", font=2);

			}

		}
	}
	mtext(page_number, side=1, outer=T, cex=.8);
}


#------------------------------------------------------------------------------

met_flux_idx_mat=find_exchange_into_shared(exch_rxn_idx, model, "RFe", "TDe")

# Print out table
if(0){
	met_ids=met_id(model);
	rxn_ids=react_name(model);
	for(i in 1:nrow(met_flux_idx_mat)){
		cat(met_ids[met_flux_idx_mat[i,1]], ",'", 
			rxn_ids[met_flux_idx_mat[i,2]], "','", 
			rxn_ids[met_flux_idx_mat[i,3]], "','",
			rxn_ids[met_flux_idx_mat[i,4]], 
			"'\n", sep="");
	}
}

# Squash matrix into vector
rxn_ids_for_fva=sort(unique(as.vector(met_flux_idx_mat)));
rxn_ids_for_fva=rxn_ids_for_fva[rxn_ids_for_fva!=0];

# Allocate space for fva results
num_indiv_exch_rxns=length(rxn_ids_for_fva);
indiv_lb_mat=matrix(0, nrow=num_indiv_exch_rxns, ncol=num_increments);
indiv_ub_mat=matrix(0, nrow=num_indiv_exch_rxns, ncol=num_increments);

for(i in 1:num_increments){

	# Prepare new stoichiometry matrix with weighted biomass proportions
	p=increments[i];
	cat("Working on ", sprintf("%03.02g", p) , " ", NameA, " vs. ", sprintf("%03.02g",(1-p)), " ", NameB, "\n", sep="");
	combined_stoich=A_stoich*p + B_stoich*(1-p);
	S[,biomassB_idx]=combined_stoich;
	S(model)=S;
	
	# Run and store FVA
	fva_out=fluxVar(model, rxn_ids_for_fva);
	indiv_lb_mat[,i]=attr(fva_out, "lp_obj")[1:num_indiv_exch_rxns];
	indiv_ub_mat[,i]=attr(fva_out, "lp_obj")[(num_indiv_exch_rxns+1):(2*num_indiv_exch_rxns)];
}


pdf(paste(Output, ".indiv.fva.pdf", sep=""),  height=11, width=8.5);
par(mfrow=c(4,1)); 
par(oma=c(1,0,1,0));
plot_flux(combined_optimal_flux, NameA, NameB, Output);
plot_individual_exch_rxns(met_flux_idx_mat, rxn_ids_for_fva, indiv_lb_mat, indiv_ub_mat, model, NameA, NameB);
dev.off();

###############################################################################

print(warnings());
cat("Done.\n");

