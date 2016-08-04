

###############################################################################

seqinr_to_matrix=function(seqinr_rec){
# Converts seqinr record into a alignment matrix
	nrow=seqinr_rec$nb
	ncol=length(strsplit(seqinr_rec$seq[[1]], "")[[1]]);	

	seq_mat=matrix("", nrow=nrow, ncol=ncol);

	for(i in 1:nrow){
		uppered=toupper(seqinr_rec$seq[[i]]);
		seq_mat[i,]=strsplit(uppered, "")[[1]];
	}

	rownames(seq_mat)=seqinr_rec$nam;
	return(seq_mat);
}

###############################################################################

clean_leadtrail_gaps=function(alignment_matrix){
# Removes leading and trailing hypens so they aren't treated as gaps
	num_samples=nrow(alignment_matrix);
	alignment_len=ncol(alignment_matrix);

	for(i in 1:num_samples){

		for(j in 1:alignment_len){
			if(alignment_matrix[i,j]=="-"){
				alignment_matrix[i,j]=" ";
			}else{
				break;
			}
		}

		for(j in alignment_len:1){
			if(alignment_matrix[i,j]=="-"){
				alignment_matrix[i,j]=" ";
			}else{
				break;
			}
		}

	}
	return(alignment_matrix);
}

###############################################################################

build_contingency_table=function(alignment_matrix, group1, group2, position){

	matrix_names=rownames(alignment_matrix);
	
	# Confirm that we have sequences in alignment
	in_g1_notinMatrix=setdiff(group1, matrix_names);
	if(length(in_g1_notinMatrix)>0){
		cat("\n\nError:  Missing names (Group 1) from alignment.\n");	
		print(in_g1_notinMatrix);
		cat("\n\n");
		quit(status=-1);
	}
	in_g2_notinMatrix=setdiff(group2, matrix_names);
	if(length(in_g2_notinMatrix)>0){
		cat("\n\nError:  Missing names (Group 2) from alignment.\n");	
		print(in_g2_notinMatrix);
		cat("\n\n");
		quit(status=-1);
	}

	align_len=ncol(alignment_matrix);
	if(position>align_len){
		cat("Error: Requested position greater than alignment length.\n");
	}

	group1_residues=alignment_matrix[group1, position];
	group2_residues=alignment_matrix[group2, position];

	# Identify positions with spaces
	g1_spaces=which(group1_residues==" ");
	g2_spaces=which(group2_residues==" ");

	# Remove spaces
	if(length(g1_spaces)){
		group1_residues=group1_residues[-g1_spaces];
	}
	if(length(g2_spaces)){
		group2_residues=group2_residues[-g2_spaces];
	}

	# Compute the residue counts for each group
	g1_counts=table(group1_residues);
	g2_counts=table(group2_residues);

	residues=union(names(g1_counts), names(g2_counts));
	num_residues=length(residues);

	if(num_residues==1){
		# If no difference in residues, return -1
		return(matrix());
	}else{
		# Create a empty contingency table
		cont_table=matrix(0, nrow=num_residues, ncol=2);
		rownames(cont_table)=residues;
		colnames(cont_table)=c("G1", "G2");

		# Populate contigency table with counts
		for(r in names(g1_counts)){
			cont_table[r,1]=g1_counts[r]
		}
		for(r in names(g2_counts)){
			cont_table[r,2]=g2_counts[r]
		}

		return(cont_table);
	}
}

###############################################################################

summarize_contingency_table=function(cont_table){

	g1max=max(cont_table[,1]);
	g1tot=sum(cont_table[,1]);

	g2max=max(cont_table[,2]);
	g2tot=sum(cont_table[,2]);

	g1maxidx=min(which(cont_table[,1]==g1max));
	g2maxidx=min(which(cont_table[,2]==g2max));
	
	residues=rownames(cont_table);
	g1_nuc=residues[g1maxidx];
	g2_nuc=residues[g2maxidx];

	g1_prop=100*g1max/g1tot;
	g2_prop=100*g2max/g2tot;
	
	return(c(
		sprintf("%5.1f", g1_prop),
		g1_nuc,
		sprintf("%5.1f", g2_prop),
		g2_nuc
		)
	);
}	

###############################################################################

load_accession_to_strainname_map=function(map_file){
# Loads access to strain names and places it in a list
	map=list();
	cat("Loading Accession-to-StrainName file: ", map_file, "\n", sep="");
	tab=as.matrix(read.table(map_file, sep="\t"));
	num_map_values=nrow(tab);
	cat("Num entries in map file: ", num_map_values, "\n");
	for(i in 1:num_map_values){
		map[[tab[i,1]]]=tab[i,2];
	}	
	return(map);
}

###############################################################################

if(0){
	map=load_accession_to_strainname_map("H1N1_2012.map");
	print(map);
}

if(0){
	library(seqinr);

	seqinr_rec=read.alignment("peptides_wReference.aln", format="clustal");

	mat=seqinr_to_matrix(seqinr_rec);
	mat=clean_leadtrail_gaps(mat);

	num_samples=nrow(mat);
	cat("Num samples: ", num_samples, "\n");

	#189
	for(i in 1:num_samples){
		ct=build_contingency_table(mat, 1:100, 150:250, i);
		if(ncol(ct)>1){
			ft=fisher.test(ct);
			pval=ft$p.val;
			if(pval<.05){
				cat(i, ": ", pval, "\n");
			}
			ct_string=summarize_contingency_table(ct);
			print(ct_string);
				
		}
	}
}




