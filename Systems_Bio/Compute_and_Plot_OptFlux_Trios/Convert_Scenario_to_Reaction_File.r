#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
		"input", "i", 1, "character",
		"scenario_col", "c", 2, "numeric"
);

DEF_SCEN_COL=10;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"\n",
		"	-i <input scenario file name>\n",
		"	[-c <column to find first scenario, def=", DEF_SCEN_COL, ">]\n",
		"\n",
		"This script will read in a scenario file and\n",
		"generate a reaction file, for each scenario.\n",
		"\n",
		"The scenario file is just a reaction file\n",
		"but with the farthest right columns filled with\n",
		"lowerbound values.\n",
		"\n",
		"The name of the scenario is in the header line (first row).\n",
		"If the row value in the scenario column is filled\n",
		"in, then that will be the value used in the scenario.\n",
		"If it is blank, then the default will be used from\n",
		"the 'lowbnd' column for that respective row.\n",
		"\n");

if(!length(opt$input)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input;
ScenarioCol=DEF_SCEN_COL;

if(length(opt$scenario_col)){
	ScenarioCol=opt$scenario_col;
}

OutputFileRoot=gsub("\\.tsv$", "", InputFile);

cat("Input File: ", InputFile, "\n", sep="");
cat("Output Root: ", OutputFileRoot, "\n", sep="");
cat("Scenario Column: ", ScenarioCol, "\n", sep="");


###############################################################################

reaction_columns=c(
	"abbreviation",
	"name",
	"equation",
	"reversible",
	"compartment",
	"lowbnd",
	"uppbnd"
);

data_matrix=as.matrix(read.table(InputFile, header=T, sep="\t"));

num_col=ncol(data_matrix);
scenarios=data_matrix[,ScenarioCol:num_col];
num_scenarios=ncol(scenarios);

cat("Num scenarios found: ", num_scenarios, "\n");

scenarios=apply(scenarios, 2, as.numeric);
scenario_names=colnames(scenarios);


apply_scenario=function(scenario, default){
	
	num_settings=length(scenario);
	num_default=length(default);
	if(num_settings != num_default){
		cat("Error: settings don't match default length.\n");
	}

	scenario_val=as.numeric(scenario);
	
	out_settings=numeric(num_settings);
	for(i in 1:num_settings){
		if(!is.na(scenario_val[i])){
			out_settings[i]=scenario_val[i];
		}else{
			out_settings[i]=default[i];
		}
	}
	
	return(out_settings);
}



for(i in 1:num_scenarios){
	
	scenario_lb=scenarios[,i];
	scenario_name=scenario_names[i];

	cat("Working on scenario: ", scenario_name, "\n");
	output_fname=paste(OutputFileRoot, ".", scenario_name, ".bnds.tsv", sep="");

	#print(data_matrix);
	output_matrix=data_matrix[,reaction_columns];

	new_lb=apply_scenario(scenario_lb, as.numeric(output_matrix[,"lowbnd"]));

	output_matrix[,"lowbnd"]=as.character(new_lb);
	
	write.table(output_matrix, file=output_fname, quote=F, sep="\t", col.names=T, row.names=F);

}






###############################################################################

print(warnings());
cat("Done.\n");


























