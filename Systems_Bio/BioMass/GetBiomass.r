
library(sybilSBML);
library(sybil);


get_biomass_met=function(model){

	# Find biomass reaction by name
	react_name=react_name(model);
	obj_func_idx=which(obj_coef(model)==1);
	if(length(obj_func_idx)!=1){
		cat("Error: Cannot identify Biomass/Objective function.\n");
		return;
	}else{
		cat("Biomass/Objective:\n");
		cat("\t", react_name[obj_func_idx], "\n", sep="");
		cat("Column:\n");
		cat("\t", obj_func_idx, "\n", sep="");
	}

	# Metabolites: Row, Reactions: Col
	S=S(model);
	biomass_metabolites=S[, obj_func_idx]
	biomass_in_met_idx=which(biomass_metabolites<0);
	biomass_out_met_idx=which(biomass_metabolites>0);

	# In Met Names
	met_names=met_name(model);

	biomass_names=list();
	biomass_names[["Inputs"]]=sort(met_names[biomass_in_met_idx]);	
	biomass_names[["Outputs"]]=sort(met_names[biomass_out_met_idx]);	
	return(biomass_names);

}

model_uploaded=readSBMLmod("/home/kli/cobra/Thiobacillus Denitrificans/RAST_Made/Seed292415.7.96972.uploaded.mod2.xml");
model_premade=readSBMLmod("/home/kli/cobra/Thiobacillus Denitrificans/TheSeed_Premade/Seed292415.3.xml");

uploaded_biomass_rxn=get_biomass_met(model_uploaded);
premade_biomass_rxn=get_biomass_met(model_premade);

