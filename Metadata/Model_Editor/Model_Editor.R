#!/usr/bin/Rscript

library(tcltk2);

# Global variables
variable_lists=list();

if(0){
	variable_lists[["removed"]]=c("barcode_id");
	variable_lists[["available"]]=c("height", "weight", "cd4");
	variable_lists[["covariates"]]=c("bmi", "sex", "race");
	variable_lists[["repeat_identifier"]]=c("date");
	variable_lists[["groups"]]=list();
	variable_lists[["groups"]][["cv"]]=c("pulse", "systolic", "diastolic");
	variable_lists[["groups"]][["pulm"]]=c("FEV1_pred", "FVC_pred", "DLCO_pred", "VA_pred");
}

list_box_list=list();
model_changed=F;
current_filename="";


load_variable_list=function(filename){
	cat("Loading Variable List: '", filename, "'\n");
	fh=file(filename, "r");

	done=F;
	variable_lists=list();
	while(!done){
		line=readLines(fh, n=1);
		if(length(line)==0){
			done=T;
		}else{
			components=strsplit(line, "\t")[[1]];

			var_list=c();
			if(length(components)==2){
				var_list=strsplit(components[2], ",")[[1]];
			}
			group_cat=strsplit(components[1], "/")[[1]];

			if(length(var_list)==0){
				var_list=c();
				cat("Empty List for Variable Category: ", components[1], "\n");
			}

			if(group_cat[1]=="groups"){
				variable_lists[["groups"]][[group_cat[2]]]=var_list;
			}else{
				variable_lists[[components[1]]]=var_list;
			}
		}
	}
	close(fh);

	# Confirm that data read in is what is expected
	listnames=names(variable_lists);
	variable_categories= c("removed", "available", "covariates", "groups", "repeat_identifier");
	excess=setdiff(listnames, variable_categories);
	expected=intersect(listnames, variable_categories);
	if(length(excess)>0 || length(expected)==0){
		tk_messageBox(type="ok", 
			caption="Loading model file...",
			message="The model file does not appear to be valid.  Try importing the TSV file headers instead?");
		return(-1);
	}

	if(is.null(variable_lists[["removed"]])){ variable_lists[["removed"]]=character();}
	if(is.null(variable_lists[["available"]])){ variable_lists[["available"]]=character();}
	if(is.null(variable_lists[["covariates"]])){ variable_lists[["covariates"]]=character();}
	if(is.null(variable_lists[["groups"]])){ variable_lists[["groups"]]=list();}
	if(is.null(variable_lists[["repeat_identifier"]])){ variable_lists[["repeat_identifier"]]=character();}

	variable_lists<<-variable_lists;
	model_changed<<-FALSE;
	current_filename<<-filename;
	cat("done.\n");
	return(0);
}

write_variable_list=function(filename){

	file_comp=strsplit(filename, "\\.")[[1]];
	if(length(file_comp)==1){
		filename=paste(filename, ".mod", sep="");
	}

	cat("Writing variable list to: ", filename, "\n");
	fh=file(filename, "w");
	
	categories=names(variable_lists);
	for(nm in categories){
		if(nm!="groups"){
			cat(file=fh, nm, "\t", paste(variable_lists[[nm]], collapse=","), "\n", sep="");
		}else{
			grp_nm=names(variable_lists[["groups"]]);
			for(gnm in grp_nm){
				cat(file=fh, paste(nm, "/", gnm, sep=""), 
				    "\t", paste(variable_lists[[nm]][[gnm]], collapse=","), "\n", sep="");
			}
		}
	}

	close(fh);
	cat("done.\n");

	model_changed<<-FALSE;
}

import_tsv_header=function(filename){
	cat("Loading TSV file for variables: '", filename, "'\n");
	fh=file(filename, "r");
	line=readLines(fh, n=1);
	close(fh);
	var_names=strsplit(line, "\t")[[1]];
	var_names=setdiff(var_names, "");
	variable_lists=list();
	variable_lists[["removed"]]=character();
	variable_lists[["available"]]=var_names;
	variable_lists[["covariates"]]=character();
	variable_lists[["groups"]]=list();
	variable_lists[["repeat_identifier"]]=character();
	variable_lists<<-variable_lists;
}

rename_group=function(old_grpname, new_grpname){
	variable_lists[["groups"]][[new_grpname]]=variable_lists[["groups"]][[old_grpname]];	
	variable_lists[["groups"]][[old_grpname]]=c();	
	variable_lists<<-variable_lists;
}


move_variables=function(src, dst, src_grp=NULL, dst_grp=NULL, var_list, win){

	cat("Moving from: ", src, " to: ", dst, "\n");

	if(dst!=""){
		if(is.null(dst_grp)){
			variable_lists[[dst]]=union(variable_lists[[dst]], var_list);
		}else{
			variable_lists[["groups"]][[dst_grp]]=union(variable_lists[["groups"]][[dst_grp]], var_list);
		}
	}else{
		cat("Destination is empty...\n");
	}

	if(is.null(src_grp)){
		variable_lists[[src]]=setdiff(variable_lists[[src]], var_list);

		# Remove from repeated identifier list too if it was a covariate
		if(src=="covariates"){
			cur=variable_lists[["repeat_identifier"]];
			upd=setdiff(cur, var_list);
			variable_lists[["repeat_identifier"]]=upd;
		}
	}else{
		cat("Source is: ", src, " (", src_grp, ") \n");
		variable_lists[["groups"]][[src_grp]]=setdiff(variable_lists[["groups"]][[src_grp]], var_list);
	}

	# Write to global variable
	variable_lists<<-variable_lists;

	if(src=="removed" || dst=="removed"){
		update_list_box("removed");
	}
	if(src=="available" || dst=="available"){
		update_list_box("available");
	}
	if(src=="covariates" || dst=="covariates"){
		update_list_box("covariates");
	}
	if(src=="groups"){
		update_list_box("groups",src_grp);
		update_list_box("available");
	}

	model_changed<<-TRUE;
	update_title(win);
}

add_to_repeated_measures=function(target_var){
	cur_ids=variable_lists[["repeat_identifier"]];
	to_add=setdiff(target_var, cur_ids);
	new_set=c(cur_ids, to_add);
	variable_lists[["repeat_identifier"]]=new_set;
	variable_lists<<-variable_lists;
}

get_variable_group_counts=function(){
	group_names=names(variable_lists[["groups"]]);
	counts_list=list();
	all=c();
	for(gn in group_names){
		all=c(all, variable_lists[["groups"]][[gn]]);
	}
	t=table(all);
	return(t);
}

remove_variable_from_groups=function(to_remove, group_name){
	grps=names(variable_lists[["groups"]]);
	for(gn in grps){
		orig_members=variable_lists[["groups"]][[gn]];
		new_members=setdiff(orig_members, to_remove);
		variable_lists[["groups"]][[gn]]=new_members;
	}

	variable_lists<<-variable_lists;
}

count_affected_groups=function(var){
	grps=names(variable_lists[["groups"]]);
	num_gps_affected=0;
	for(gn in grps){
		members=variable_lists[["groups"]][[gn]];
		overlap=intersect(members, var);
		if(length(overlap)>0){
			num_gps_affected=num_gps_affected+1;
		}
	}
	return(num_gps_affected);

}

delete_group=function(grp_name, win){

	variable_lists[["groups"]][[grp_name]]=NULL;
	variable_lists<<-variable_lists;
	x=(tkwinfo("x", win));
	y=(tkwinfo("y", win));
	h=(tkwinfo("height", win));
	w=(tkwinfo("width", win));
	loc=paste("+",x,"+",y, sep="");
	tkdestroy(win);
	model_changed<<-TRUE;

	num_grps=length(variable_lists[["groups"]]);
	new_window_width=as.integer(as.integer(w)*(3+num_grps)/(3+num_grps+1));

	win_main<<-draw_main(loc, h, new_window_width);
}

update_title=function(win){
	chg_flg=ifelse(model_changed,"*", "");
	fname=ifelse(current_filename=="", "[untitled]", current_filename);
	title=paste(fname, chg_flg, " - CMM Model Edtor");
	tktitle(win)=title;
}

draw_remove_ctls=function(win, row, col){

	onButt=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)));
		move_variables("removed", "available", NULL, NULL, variable_lists[["removed"]][select_ix+1], win);
	}

	listbox_h=tk2listbox(win, height = 10, selectmode = "extended");

	tkgrid(tk2label(win, text="Removed", justify="left"), row=1+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(listbox_h, row=2+row, column=col, 
		padx=c(15,5), pady=c(5,5), sticky="nsew");
	tkgrid(tk2button(win, text = "Make Available", width = -6, command = onButt), row=3+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	return(listbox_h);
}

add_and_select_groups_dialog=function(variables_to_add){
	
	dlg=tktoplevel();
	num_groups_added=0;

	onAddNewGroup=function(){
		avail_groups=names(variable_lists[["groups"]]);
		new_grp_name=tclvalue(newGroupName);
		if(new_grp_name==""){ return();}

		if(any(new_grp_name==avail_groups)){
			tk_messageBox(type="ok", 
				caption="Error adding new group",
				message=paste(new_grp_name, " already exists as a group."));
			
		}else{
			variable_lists[["groups"]][[new_grp_name]]=character();
			variable_lists<<-variable_lists;

			tkinsert(avail_group_listbox, "end", new_grp_name);
			tkdelete(new_group_name_entry,0,"end");
			num_groups_added<<-num_groups_added+1;
		}

	}
	onDone=function(){
		tkdestroy(dlg);
	}
	onAddToGroup=function(){
		avail_groups=names(variable_lists[["groups"]]);
		cat("Avail Groups:\n");
		print(avail_groups);
		select_ix=(as.numeric(tkcurselection(avail_group_listbox)))+1;
		if(length(select_ix)==0){
			tk_messageBox(type="ok", 
				caption="Error adding variables to group",
				message=paste("Please select a group first."));
			
		}else{
			dst_groups=avail_groups[select_ix];
			for(gn in dst_groups){
				oldlist=variable_lists[["groups"]][[gn]];
				toadd=setdiff(variables_to_add, oldlist);
				newlist=c(oldlist, toadd);
				variable_lists[["groups"]][[gn]]=newlist;
			}
			variable_lists<<-variable_lists;
			tkdestroy(dlg);
		}
	}
	
	add_new_group_label=tk2label(dlg, text="Create New Group:");
	newGroupName=tclVar("");
	new_group_name_entry=tk2entry(dlg, width="25", textvariable=newGroupName);
	add_group_button=tk2button(dlg, text="Create", width=-6, command=onAddNewGroup);
	avail_group_name_label=tk2label(dlg, text="Currently Available Groups:");
	avail_group_listbox=tk2listbox(dlg, height=10, selectmode="extended", scroll="both");
	add_to_groups_button=tk2button(dlg, text="Add Variables to Group(s)", width=-6, command=onAddToGroup);
	done_button=tk2button(dlg, text="Ok", width=-6, command=onDone);

	tkgrid(add_new_group_label, sticky="w", pady=c(10,0));
	tkgrid(new_group_name_entry, add_group_button, padx=c(10,10), sticky="ew");
	tkgrid(avail_group_name_label, sticky="w", pady=c(10,0));
	tkgrid(avail_group_listbox, add_to_groups_button, padx=c(10,10), sticky="ew");
	tkgrid(done_button, pady=c(30,20), columnspan=2);

	avail_groups=names(variable_lists[["groups"]]);
	for(gn in avail_groups){
		tkinsert(avail_group_listbox, "end", gn);
	}

	tkwm.deiconify(dlg);
	tktitle(dlg)="Add variables to groups...";
	tkwm.resizable(dlg, 0, 0);
	tkgrab.set(dlg);
	tkfocus(dlg);
	tkwait.window(dlg);

	return(num_groups_added);

}

draw_available_ctls=function(win, row, col){

	onButt_mkcov=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)))+1;
		target_var=variable_lists[["available"]][select_ix];

		num_affected=count_affected_groups(target_var);
		if(num_affected>0){
			msg=paste("Changing these variables into covariates will remove them from ", num_affected, " groups.", sep="");
			res=tkmessageBox(type="okcancel", message=msg);
			res_str=tclvalue(res);
			if(res_str!="ok"){
				return();
			}
			remove_variable_from_groups(target_var);
		}
		move_variables("available", "covariates", NULL, NULL, target_var, win);
		update_list_box("groups");
	}
	onButt_grpadd=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)))+1;
		if(length(select_ix)==0){
			tk_messageBox(type="ok",
			      message="Please select some variables first.");
			return();
		}
		# query user for group
		num_groups_added=add_and_select_groups_dialog(variable_lists[["available"]][select_ix]);

		if(num_groups_added){
			cat("Refreshin window because group added...\n");
			x=(tkwinfo("x", win));
			y=(tkwinfo("y", win));
			loc=paste("+",x,"+",y, sep="");
			h=(tkwinfo("height", win));
			w=(tkwinfo("width", win));
			tkdestroy(win);

			num_grps=length(variable_lists[["groups"]]);
			new_width=as.integer(as.integer(w)*(3+num_grps)/(3+num_grps-num_groups_added));
			win_main<<-draw_main(loc, h, new_width);
		}
		model_changed<<-TRUE;
		update_list_box("_all");
	}
	onButt_repid=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)))+1;
		target_var=variable_lists[["available"]][select_ix];

		num_affected=count_affected_groups(target_var);
		if(num_affected>0){
			msg=paste("Making these variables repeated measure identifiers will also remove them from ", num_affected, " groups.", sep="");
			res=tkmessageBox(type="okcancel", message=msg);
			res_str=tclvalue(res);
			if(res_str!="ok"){
				return();
			}
			remove_variable_from_groups(target_var);
		}
		remove_variable_from_groups(target_var);
		move_variables("available", "covariates", NULL, NULL, target_var, win);
		add_to_repeated_measures(target_var);
		update_list_box("groups");
		update_list_box("covariates");
	}
	onButt_rm=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)))+1;
		target_var=variable_lists[["available"]][select_ix];

		num_affected=count_affected_groups(target_var);
		if(num_affected>0){
			msg=paste("Removing these variables will also remove them from ", num_affected, " groups.", sep="");
			res=tkmessageBox(type="okcancel", message=msg);
			res_str=tclvalue(res);
			if(res_str!="ok"){
				return();
			}
			remove_variable_from_groups(target_var);
		}
		move_variables("available", "removed", NULL, NULL, target_var, win);
		update_list_box("groups");
	}

	onSelect=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)))+1;
		target_var=variable_lists[["available"]][select_ix];
		set_selection(target_var);
	}

	listbox_h=tk2listbox(win, height = 10, selectmode = "extended", scroll="both");

	tkbind(listbox_h, "<<ListboxSelect>>", onSelect);

	tkgrid(tk2label(win, text="Available", justify="left"), row=1+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(listbox_h, row=2+row, column=col, 
		padx=c(5,5), pady=c(5,5), sticky="nsew");

	tkgrid(tk2button(win, text = "Make Covariate", width = -6, command = onButt_mkcov), row=3+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(tk2button(win, text = "Add to Group", width = -6, command = onButt_grpadd), row=4+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(tk2button(win, text = "Set Repeat Identifier", width = -6, command = onButt_repid), row=5+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(tk2button(win, text = "Remove", width = -6, command = onButt_rm), row=6+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	return(listbox_h);
}

draw_covariates_ctls=function(win, row, col){

	onButt=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)));
		move_variables("covariates", "available", NULL, NULL, variable_lists[["covariates"]][select_ix+1], win);
	}

	listbox_h=tk2listbox(win, height = 10, selectmode = "extended");
	
	tkgrid(tk2label(win, text="Covariates", justify="left"), row=1+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	tkgrid(listbox_h, row=2+row, column=col, 
		padx=c(5,5), pady=c(5,5), sticky="nsew");
	tkgrid(tk2button(win, text = "Make Available", width = -6, command = onButt), row=3+row, column=col, 
		padx=c(5,5), pady=c(5,5));
	return(listbox_h);
}



draw_groups_ctls=function(group_name, win, row, col){

	onRemButt=function(){
		select_ix=(as.numeric(tkcurselection(listbox_h)));
		move_variables("groups", "", group_name, NULL, variable_lists[["groups"]][[group_name]][select_ix+1], win);
	}

	onDeleteGroupButt=function(){
		num_var_in_grp=length(variable_lists[["groups"]][[group_name]]);
		if(num_var_in_grp>0){
			msg=paste("Are you sure you want to delete this group?\nIt contains ", num_var_in_grp, " variable(s)", sep="");
			res=tkmessageBox(type="okcancel", message=msg);
			res_str=tclvalue(res);
			if(res_str!="ok"){
				return();
			}
		}
		delete_group(group_name, win);
		update_list_box("_all");
	}

	onClick=function(){

		onFocusOut=function(){
			new_name=tclvalue(gname);
			cat("New Name: ", new_name, "\n");

			cur_grps=names(variable_lists[["groups"]]);
			if(new_name==group_name){
				cat("No change to group name.\n");
			}
			else if(any(new_name==cur_grps)){
				msg=paste("This group name already exists.  Reverting to original name.");
				res=tkmessageBox(type="ok", message=msg);
			}else{
				rename_group(group_name, new_name);
				list_box_list[["groups"]][[new_name]]=list_box_list[["groups"]][[group_name]];
				list_box_list[["groups"]][[group_name]]=c();	
				list_box_list<<-list_box_list;
				model_changed<<-T;
				update_title(win);
				tkconfigure(grp_label_h, text=paste("Group: ", new_name, sep=""));
			}
			tkdestroy(entry_h);
		}

		cat("Group name clicked: ", group_name, "\n");
		gname=tclVar(group_name);
		entry_h=tk2entry(win, textvariable=gname);

		tkfocus(entry_h);
		tkgrid(entry_h, row=1+row, column=col);
		tkbind(entry_h, "<FocusOut>", onFocusOut);
		tkbind(entry_h, "<Return>", onFocusOut);
	}


	listbox_h=tk2listbox(win, height = 10, selectmode = "extended");

	grp_label_h=tk2label(win, text=paste("Group: ", group_name, sep=""), justify="left");
	tkbind(grp_label_h, "<Button-1>", onClick);

	tkgrid(
	       	grp_label_h,
		row=1+row, column=col,
		padx=c(5,5), pady=c(5,5)
		);
	tkgrid(listbox_h, row=2+row, column=col, 
		padx=c(5,5), pady=c(5,5), sticky="nsew");
	tkgrid(tk2button(win, text = "Remove from Group", width = -6, command = onRemButt), row=3+row, column=col, 
		padx=c(5,5), pady=c(5,5));

	tkgrid(tk2button(win, text = "Delete Group", width = -6, command = onDeleteGroupButt), row=5+row, column=col, 
		padx=c(5,5), pady=c(5,5));


	return(listbox_h);
}

save_changes_check=function(reason){
		res=tkmessageBox(message=paste("Save changes before ", reason, "?", sep=""),
				 icon="question", type="yesnocancel", default="yes"); 
		res_str=tclvalue(res);
		if(res_str=="no"){
			return("ok");
		}else if(res_str=="yes"){
			if(current_filename==""){
				filename=tclvalue(tkgetSaveFile());
			}
			write_variable_list(current_filename);
			return("ok");
		}else if(res_str=="cancel"){
			return("cancel");
		}
}


init_menus=function(win){
	onOpenExistingModel=function(){

		if(model_changed){
			res=save_changes_check("opening another model");
			if(res=="cancel"){
				return();
			}
		}

		filename=tclvalue(tkgetOpenFile(filetypes="{{Model Files} {.mod}}"));
		if(filename!=""){ # cancel/noop 
			status=load_variable_list(filename);
			if(status==0){
				x=(tkwinfo("x", win));
				y=(tkwinfo("y", win));
				loc=paste("+",x,"+",y, sep="");
				h=(tkwinfo("height", win));
				w=(tkwinfo("width", win));
				tkdestroy(win);
				model_changed<<-FALSE;
				win_main<<-draw_main(loc, h, w);
				current_filename<<-filename;
				update_list_box("_all");
			}
		}
		cat("Open Model File Selection: ", filename, "\n");
	}
	onOpenTSVMetadataFile=function(){

		if(model_changed){
			res=save_changes_check("importing new variable list");
			if(res=="cancel"){
				return();
			}
		}

		filename=tclvalue(tkgetOpenFile(filetypes="{{Metadata TSV Files} {.tsv .txt}}"));
		if(filename!=""){ # cancel/noop 
			import_tsv_header(filename);
			x=(tkwinfo("x", win));
			y=(tkwinfo("y", win));
			loc=paste("+",x,"+",y, sep="");
			h=(tkwinfo("height", win));
			w=(tkwinfo("width", win));
			tkdestroy(win);
			current_filename<<-"";
			model_changed<<-TRUE;
			win_main<<-draw_main(loc, h, w);
			update_list_box("_all");
		}
		cat("Open TSV File Selection: ", filename, "\n");
	}
	onSaveModel=function(){
		if(current_filename!=""){ 
			write_variable_list(current_filename);
		}else{
			filename=tclvalue(tkgetSaveFile(initialfile="model", filetypes="{{Model Files} {.mod}}"));
			write_variable_list(filename);
			current_filename<<-filename;
		}
		model_changed<<-FALSE;
		update_title(win);
	}
	onSaveModelAs=function(){
		filename=tclvalue(tkgetSaveFile(initialfile="model", filetypes="{{Model Files} {.mod}}"));
		if(filename==""){ # cancel/noop 
			return();
		}
		cat("Save Model As Selection: ", filename, "\n");
		write_variable_list(filename);
		current_filename<<-filename;
		model_changed<<-FALSE;
		update_title(win);
	}
	onQuit=function(){
		if(model_changed){
			res=save_changes_check("quiting");
			if(res=="cancel"){
				return();
			}
		}
		tkdestroy(win);
	}
	onHelp=function(){
		helplist=c(
			"Choose variables from the 'Available' column.",
			"",
			"Available variables may be removed from the model, or may be assigned as 'Covariates', or assigned to multiple 'Groups'.",
			"",
			"Covariates are kept as predictors in all models, e.g. sex, age, race, etc.",
			"Variables moved into 'Groups' will be analyzed together and may be used as predictors or multivariate responses.",
			"'Removed' variables are just set aside from a particular analysis.",
			"",
		       	"An 'Available' variable may be assigned to multiple groups, but a 'Covariate' and 'Removed' variable may not be.",
			"",
			"Variables may be assigned to different categories by selecting them from the list and clicking on the appropriate buttons below the list that contain them.",
			"",
			"Setting a variable as a repeat(ed measure) identifier removes the variable from being available.",
			"Repeated measure identifiers are usually identifier variables that differentiate samples from the same subject, e.g. time",
			"",
			"Remember to save your model!",
			"You can start a new model by importing the variables names from your TSV formatted metadata file.",
			"Previously saved models may also be edited and viewed by 'Opening Existing Model...'"
		);
		tkmessageBox(message=paste(helplist, collapse="\n"), type="ok");
	}

	onAbout=function(){
		tkmessageBox(message="CMM Model Editor\nVersion 0.5 (beta)\nCenter for Medicine and the Microbiome\nCopyright 2018", type="ok");
	}

	menu_obj=tk2menu(win);
	tkconfigure(win, menu=menu_obj);
	file_menu=tk2menu(menu_obj, tearoff=FALSE);
	tkadd(file_menu, "command", label="Open Existing Model...", command=onOpenExistingModel);
	tkadd(file_menu, "command", label="Import TSV Metadata File Headers...", command=onOpenTSVMetadataFile);
	tkadd(file_menu, "separator");
	tkadd(file_menu, "command", label="Save Model", command=onSaveModel);
	tkadd(file_menu, "command", label="Save Model As...", command=onSaveModelAs);
	tkadd(file_menu, "separator");
	tkadd(file_menu, "command", label="Quit", command=onQuit);
	tkadd(menu_obj, "cascade", label="File", menu=file_menu);


	help_menu=tk2menu(menu_obj, tearoff=FALSE);
	tkadd(help_menu, "command", label="Help", command=onHelp);
	tkadd(help_menu, "separator");
	tkadd(help_menu, "command", label="About", command=onAbout);
	tkadd(menu_obj, "cascade", label="Help", menu=help_menu);
}

set_selection=function(target_var, which_textbox){

	#cat("Selected variables:\n");
	#print(target_var);

	#cat("Available Groups:\n");
	group_names=names(list_box_list[["groups"]]);

	for(grp in group_names){
		var_list=variable_lists[["groups"]][[grp]];
		tkselection.clear(list_box_list[["groups"]][[grp]], 0, length(var_list)-1);
		for(tar in target_var){
			idx=which(tar==var_list);
			if(length(idx)){
				tkselection.set(list_box_list[["groups"]][[grp]], idx-1);
			}
		}
	}
}

update_list_box=function(target, group_name=NULL){
	if(target=="_all"){
		update_list_box("removed");
		update_list_box("available");
		update_list_box("covariates");
		update_list_box("groups");
	}else if(target=="groups" && is.null(group_name)){
		for(g in names(variable_lists[["groups"]])){
			update_list_box("groups", g);
		}
	}else if(target!="groups"){
		tkdelete(list_box_list[[target]],0, "end");
		if(target=="available"){
			count_table=get_variable_group_counts();
			avail_list=variable_lists[["available"]];
			num_avail=length(avail_list);
			if(num_avail>0){
				for(i in 1:num_avail){
					var=avail_list[i];
					count=count_table[var];
					var_wcount=ifelse(is.na(count), 
						"     ", paste("[",count,"]",sep=""));
					var_wcount=paste(var_wcount, "  ", var, sep="");
					tkinsert(list_box_list[["available"]], "end", var_wcount);
				}
			}
		}else if(target=="covariates"){
			cov_list=variable_lists[["covariates"]];
			rep_id_list=variable_lists[["repeat_identifier"]];
			num_cov=length(cov_list);
			if(num_cov>0){
				for(i in 1:num_cov){
					var=cov_list[i];
					is_rpm=any(var==rep_id_list);
					var_rp_tag=ifelse(is_rpm, "[R]", "     ");
					var_wtag=paste(var_rp_tag, "  ", var, sep="");
					tkinsert(list_box_list[["covariates"]], "end", var_wtag);
				}	
			}	
		}else{
			for(val in variable_lists[[target]]){
				tkinsert(list_box_list[[target]], "end", val);
			}
		}
	}else{
		tkdelete(list_box_list[["groups"]][[group_name]],0, "end");
		for(val in variable_lists[["groups"]][[group_name]]){
			tkinsert(list_box_list[["groups"]][[group_name]], "end", val);
		}
	}
}
	

draw_main=function(location=NULL, height=NULL, width=NULL){

	win_main<<-tktoplevel();

	if(!is.null(height) && !is.null(width)){
		tkwm.geometry(win_main, paste(width, "x", height, sep=""));
	}

	update_title(win_main);

	tkwm.resizable(win_main, 1, 1);
	if(!is.null(location)){
		tkwm.geometry(win_main, location);
	}

	init_menus(win_main);

	tkgrid(tk2label(win_main, text="CMM Model Editor", font=c("bold", 24), justify="left"), row=0, column=0,
                padx=c(5,5), pady=c(5,5), columnspan=3);

	list_box_list[["removed"]]=draw_remove_ctls(win_main, row=1, col=0);
	#readline("Press enter to continue.");
	list_box_list[["available"]]=draw_available_ctls(win_main, row=1, col=1);
	#readline("Press enter to continue.");
	list_box_list[["covariates"]]=draw_covariates_ctls(win_main, row=1, col=2);

	cat("Drawing group controls...\n");
	group_names=names(variable_lists[["groups"]]);
	grp_id=0;
	for(grp in group_names){
		list_box_list[["groups"]][[grp]]=
			draw_groups_ctls(grp, win_main, row=1, col=3+grp_id);
		grp_id=grp_id+1;
	}

	tkgrid(tk2label(win_main, text=""), row=9, column=0,
                padx=c(5,5), pady=c(5,5), columnspan=3);

	# Allow list boxes to follow size of window dimensions
	tkgrid.rowconfigure(win_main, 3, weight=1);
	num_groups=length(group_names);
	for(i in 1:(3+num_groups)){
		tkgrid.columnconfigure(win_main, i-1, weight=1);
	}

	# Set min height and width based on original spacing with adjustments to number of groups
	if(is.null(height) && is.null(width)){
		glob_min_height<<-as.numeric(tkwinfo("height", win_main));
		glob_min_width<<-as.numeric(tkwinfo("width", win_main));
	}
	min_width_wgroups=as.integer((glob_min_width/3.0)*(3+num_groups));
	tkwm.minsize(win_main, min_width_wgroups, glob_min_height);


	list_box_list<<-list_box_list;

	cat("Done with draw_main...\n");
	
}


draw_main();





