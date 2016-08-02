#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('seqinr');

params=c(
	"pathway_table", "t", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-t <pathway table>\n",
	"\n",
        "\n", sep="");

if(!length(opt$pathway_table)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------

PathwayTable=opt$pathway_table;
cat("Pathway Table Filename: ", PathwayTable, "\n");

#----------------------------------------------------------

# Load file into table
table=as.matrix(read.table(PathwayTable, sep="\t", quote="", header=T));
nrows=nrow(table);
cat("Num rows loaded: ", nrows, "\n");

# Convert numeric entries in table into integers
ids=matrix(0, nrow=nrows, ncol=3);
ids[,1]=as.numeric(table[,"id"]);
ids[,2]=as.numeric(table[,"parent_id"]);
ids[,3]=as.numeric(table[,"child_count"]);
colnames(ids)=c("id", "parent_id", "child_count");

# Identify which id's have EC, ie. leave nodes
leaves_wEC=table[,"ec_id"]!="";
num_leaves=sum(leaves_wEC);
cat("Num ECs (leaves): ", num_leaves, "\n");
leaves_wEC_ix=which(table[,"ec_id"]!="");

# Generate mapping from EC string to its multiple id's
EC_mapping=list();
for(i in leaves_wEC_ix){
	if(is.null(EC_mapping[[table[i,"ec_id"]]])){
		EC_mapping[[table[i,"ec_id"]]]=ids[i,"id"];		
	}else{
		EC_mapping[[table[i,"ec_id"]]]=c(
			EC_mapping[[table[i,"ec_id"]]], ids[i,"id"]);		
	}
}
#print(EC_mapping);

# Generate mapping from child to parent
child_parent_mapping=list();
id_name_mapping=list();

id_level_mapping=list();
id_level_backup_mapping=list();

for(i in 1:nrows){
	child_parent_mapping[[ids[i,"id"]]]=ids[i,"parent_id"];
	id_name_mapping[[ids[i,"id"]]]=table[i,"name"];

	id_level_mapping[[ids[i,"id"]]]=table[i,"level"];
	id_level_backup_mapping[[ids[i,"id"]]]=table[i,"level_backup"];
}
#print(child_parent_mapping);
#print(id_name_mapping);

travel_up_hierarchy=function(target_ec, EC_mapping, child_parent_mapping, id_name_mapping){
	# Get ID of EC
	ec_ids=EC_mapping[[target_ec]];
	cat("EC: ", target_ec, ": ", ec_ids, "\n");
	
	for(ec_id in ec_ids){
		cat("EC: ", ec_id, "\n");
		cur_id=ec_id;
		cat("\t");
		while(cur_id!=1){
			#cat("\t", id_name_mapping[[cur_id]], "\n", sep="");
			#cat("\t  ", id_level_mapping[[cur_id]], "\n", sep="");
			#cat("\t  ", id_level_backup_mapping[[cur_id]], "\n", sep="");

			cat(id_name_mapping[[cur_id]], "(", id_level_mapping[[cur_id]], ")::", sep="");
			cur_id=child_parent_mapping[[cur_id]];
		}
		cat("\n");
	}
}

target_ec="2.2.1.2";
travel_up_hierarchy(target_ec, EC_mapping, child_parent_mapping, id_name_mapping);



cat("Done.\n");
