
library(shiny)

ModelBuilder.GroupControl=function(id){
	
	tagList(
	
		textInput(
			paste("ModelBuilderTab.group.name_", id, sep=""),
			label=paste("Group ", id, ": ", sep=""),
			paste("MyGroup", id, sep="")),

		selectInput(
			paste("ModelBuilderTab.group.selector_", id, sep=""),
			NULL,
			choices="", multiple=T, selectize=F, size=12),
		
		actionButton(
			paste("ModelBuilderTab.group.delete_", id, sep=""),
			"Delete Selected", icon=icon("times")),
			
		actionButton(
			paste("ModelBuilderTab.group.clear_", id, sep=""),
			"Clear All Variables", icon=icon("trash-alt")),
			
		actionButton(
			paste("ModelBuilderTab.group.copyto_", id, sep=""),
			"Copy Selected To...", icon=icon("copy"))

	)
}

ModelBuilderTab.init_variables=function(){
	model_data=list();
	model_data[["ModelBuilderTab.excluded_selector"]]=character(0);
	model_data[["ModelBuilderTab.available_selector"]]=character(0);
	model_data[["ModelBuilderTab.covariates_selector"]]=character(0);
	return(model_data);
}

ModelBuilderTab.copy_to_dialog=function(input, session, source_grp_id){
	
	group_names=c();
	for(gid in 1:ModelBuilderTab.num_groups){
		if(gid!=source_grp_id){
			group_names=c(group_names,input[[paste("ModelBuilderTab.group.name_", gid, sep="")]]);
		}
	}

	showModal(
		modalDialog(
			selectInput("ModelBuilderTab.copy_to_dialog.group_select",
				label="Select group(s) to copy variables to:",
				choices=group_names,
				selected="",
				size=6,
				multiple=T,
				selectize=F
			),
			footer=fluidRow(
				column(4),
				column(3, actionButton("ModelBuilderTab.copy_to_dialog.okButton", label="Copy")),
				column(3, actionButton("ModelBuilderTab.copy_to_dialog.cancelButton", label="Cancel"))
				),
			size="s"
		)
	);
	
}

ModelBuilderTab=function(){

	selector_display_length=12;
	ModelBuilderTab.num_groups=6;
	grp_names=paste("untitled", 1:ModelBuilderTab.num_groups);

	tabPanel("Model Builder",
		fluidRow(
			column(3, tags$h1("Model Builder"), align="right"),
			tags$br(),
			actionLink("ModelBuilderTab.MainHelp", NULL, icon=icon("question-circle"))
		),
		tags$hr(),
		fixedRow(
			column(2,
				tags$b("Excluded"),
				actionLink("ModelBuilderTab.ExcludedHelp", NULL, icon=icon("question-circle")),
				selectInput("ModelBuilderTab.excluded_selector", 
					NULL, choices="", multiple=T, selectize=F, size=selector_display_length)
			),
			column(1, 
					tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
					actionButton("ModelBuilderTab.exclude_button",NULL, icon=icon("arrow-left")),
					actionButton("ModelBuilderTab.include_button",NULL, icon=icon("arrow-right"))		
			),
			column(3,
				tags$b("Available"),
				actionLink("ModelBuilderTab.AvailableHelp", NULL, icon=icon("question-circle")),
				selectInput("ModelBuilderTab.available_selector", 
					NULL, choices=ModelBuilder.available_variables, multiple=T, selectize=F, size=selector_display_length-4),
				tags$h2(""),
				tags$b("Select group(s) to add to:"),
				actionLink("ModelBuilderTab.GroupHelp", NULL, icon=icon("question-circle")),
				fixedRow(
					column(8,
						selectInput("ModelBuilderTab.group_selector", 
							NULL, choices=grp_names, multiple=T, selectize=T)),
					column(1,
						actionButton("ModelBuilderTab.add_to_group", "Add", icon=icon("arrow-down")))
				)
			),
			column(1,
					tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
					actionButton("ModelBuilderTab.remove_from_covar_button", NULL, icon=icon("arrow-left")),
					actionButton("ModelBuilderTab.add_to_covar_button", NULL, icon=icon("arrow-right"))		
			),
			column(2,
				tags$b("Covariates"),
				actionLink("ModelBuilderTab.CovariatesHelp", NULL, icon=icon("question-circle")),
				selectInput("ModelBuilderTab.covariates_selector", 
					NULL, choices="", multiple=T, selectize=F, size=selector_display_length)
			)
		),
		tags$hr(),
		fluidRow(
			column(2, ModelBuilder.GroupControl(1)),
			column(2, ModelBuilder.GroupControl(2)),
			column(2, ModelBuilder.GroupControl(3)),
			column(2, ModelBuilder.GroupControl(4)),
			column(2, ModelBuilder.GroupControl(5)),
			column(2, ModelBuilder.GroupControl(6))
		),
			
		uiOutput("StudyTab.group_controls")
	)

}

#------------------------------------------------------------------------------

ModelBuilderTab.move_selected=function(input, session, from, to){
		
	variables_selected=input[[from]];
	
	session$userData[["ModelBuilderSets"]][[from]]=
		sort(setdiff(session$userData[["ModelBuilderSets"]][[from]], variables_selected));
		
	session$userData[["ModelBuilderSets"]][[to]]=
		sort(union(session$userData[["ModelBuilderSets"]][[to]], variables_selected));
	
	updateSelectInput(session, from, choices=session$userData[["ModelBuilderSets"]][[from]]);
	updateSelectInput(session, to, choices=session$userData[["ModelBuilderSets"]][[to]]);
	
	cat("Leaving move_selected.\n");
}

#------------------------------------------------------------------------------
# Group Selector Functions

ModelBuilderTab.get_all_group_names=function(input){
	grp_names=c();
	for(i in 1:ModelBuilderTab.num_groups){
		control_ids=paste("ModelBuilderTab.group.name_", i, sep="");
		grp_names=c(grp_names, input[[control_ids]]);
	}
	return(grp_names);
}

ModelBuilderTab.update_grp_selector=function(input, session){
	grp_names=ModelBuilderTab.get_all_group_names(input);
	updateSelectizeInput(session, "ModelBuilderTab.group_selector", choices=grp_names);
}

ModelBuilderTab.delete_grp_var=function(input, session, id){
	ctl_name=paste("ModelBuilderTab.group.selector_", id, sep="");
	variables_selected=input[[ctl_name]];
	
	session$userData[["ModelBuilderSets"]][[ctl_name]]=
		sort(setdiff(session$userData[["ModelBuilderSets"]][[ctl_name]], variables_selected));
	updateSelectInput(session, ctl_name, choices=session$userData[["ModelBuilderSets"]][[ctl_name]]);

	#ModelBuilderTab.sets[[ctl_name]]<<-sort(setdiff(ModelBuilderTab.sets[[ctl_name]], variables_selected));
	#updateSelectInput(session, ctl_name, choices=ModelBuilderTab.sets[[ctl_name]]);
}

ModelBuilderTab.clear_grp_var=function(input, session, id){
cat("clear\n");
	ctl_name=paste("ModelBuilderTab.group.selector_", id, sep="");
	variables_selected=input[[ctl_name]];
	session$userData[["ModelBuilderSets"]][[ctl_name]]=character(0);
	updateSelectInput(session, ctl_name, choices=session$userData[["ModelBuilderSets"]][[ctl_name]]);
}

ModelBuilderTab.get_all_grouped_variables=function(session){
	
	all_grp_var=c();
	for(gix in 1:ModelBuilderTab.num_groups){
		selector_name=paste("ModelBuilderTab.group.selector_", gix, sep="");
		all_grp_var=c(all_grp_var, session$userData[["ModelBuilderSets"]][[selector_name]]);
	}
	
	return(all_grp_var);
}

ModelBuilderTab.remove_variables_from_all_grouped_variables=function(session, var_to_remove){
	
	cat("inside remove_variables_from_all_grouped_variables.\n");
	for(gix in 1:ModelBuilderTab.num_groups){
		selector_name=paste("ModelBuilderTab.group.selector_", gix, sep="");
		
		session$userData[["ModelBuilderSets"]][[selector_name]]=
			setdiff(
				session$userData[["ModelBuilderSets"]][[selector_name]],
				var_to_remove
				);
		updateSelectInput(
			session, 
			selector_name, 
			choices=session$userData[["ModelBuilderSets"]][[selector_name]]
		);
	}
	
}

#------------------------------------------------------------------------------

observe_ModelBuilderTabEvents=function(input, output, session){

	observeEvent(input$ModelBuilderTab.exclude_button,{
	
		variables_selected=input[["ModelBuilderTab.available_selector"]];
		
		grouped_variables=ModelBuilderTab.get_all_grouped_variables(session);
		in_group_var=intersect(grouped_variables, variables_selected)

		if(length(in_group_var)){
			showModal(
				modalDialog(
					"You are moving some variables already assigned to a group to the excluded list.",
					"To do this, each will be removed from all the groups they belong to.  Is this ok?",
					title="Removing variables from groups",
					footer=fluidRow(
						column(4),
						column(3, actionButton("ModelBuilderTab.remove_from_available_cancelButton", label="Cancel")),
						column(3, actionButton("ModelBuilderTab.remove_from_available_okButton", label="Ok"))
					),
				size="m"
				)
			)
		}else{	
			ModelBuilderTab.move_selected(input, session,
				"ModelBuilderTab.available_selector",
				"ModelBuilderTab.excluded_selector");
		}
	
	});
	
	observeEvent(input$ModelBuilderTab.include_button,{
		ModelBuilderTab.move_selected(input, session,
			"ModelBuilderTab.excluded_selector",
			"ModelBuilderTab.available_selector");
	});
	
	observeEvent(input$ModelBuilderTab.add_to_covar_button,{
		cat("From Covariates button pushed.\n");
		
		variables_selected=input[["ModelBuilderTab.available_selector"]];
		
		grouped_variables=ModelBuilderTab.get_all_grouped_variables(session);
		in_group_var=intersect(grouped_variables, variables_selected)

		if(length(in_group_var)){
			showModal(
				modalDialog(
					"You are moving some variables already assigned to a group to the covariates list.",
					"To do this, each will be removed from all the groups they belong to.  Is this ok?",
					title="Removing variables from groups",
					footer=fluidRow(
						column(4),
						column(3, actionButton("ModelBuilderTab.remove_from_available_cancelButton", label="Cancel")),
						column(3, actionButton("ModelBuilderTab.remove_from_available_okButton", label="Ok"))
					),
				size="m"
				)
			)
		}else{	
			ModelBuilderTab.move_selected(input, session,
				"ModelBuilderTab.available_selector",
				"ModelBuilderTab.covariates_selector");
		}
		
		cat("Leaving observeEvent add to covar button\n");
	});
	
	observeEvent(input$ModelBuilderTab.remove_from_available_cancelButton, {
		removeModal();
	});
	
	observeEvent(input$ModelBuilderTab.remove_from_available_okButton, {
	
		ModelBuilderTab.remove_variables_from_all_grouped_variables(
			session, input[["ModelBuilderTab.available_selector"]]);
		
		removeModal();
			
		if(!is.null(input$ModelBuilderTab.add_to_covar_button)){
			ModelBuilderTab.move_selected(input, session,
				"ModelBuilderTab.available_selector",
				"ModelBuilderTab.covariates_selector");
		}
		
		if(!is.null(input$ModelBuilderTab.exclude_button)){
			ModelBuilderTab.move_selected(input, session,
				"ModelBuilderTab.available_selector",
				"ModelBuilderTab.excluded_selector");
		}
		
	});
	
	
	observeEvent(input$ModelBuilderTab.remove_from_covar_button,{
		ModelBuilderTab.move_selected(input, session,
			"ModelBuilderTab.covariates_selector",
			"ModelBuilderTab.available_selector");
	});

	#--------------------------------------------------------------------------

	observeEvent(input$ModelBuilderTab.add_to_group, {
		cat("Add to group button pushed.\n");
		
		all_grp_names=ModelBuilderTab.get_all_group_names(input);
		#cat("All current group names:\n");
		#print(all_grp_names);
		
		grps_selected=input$ModelBuilderTab.group_selector;
		num_grps_selected=length(grps_selected);
		#cat("Selected Groups:\n");
		#print(grps_selected);
		
		variables_selected=input$ModelBuilderTab.available_selector;
		num_var_selected=length(variables_selected)
		#cat("Selected Variables:\n");
		#print(variables_selected);
		
		# Add to group
		selected_grps=c();
		if(num_grps_selected==0 || num_var_selected==0){
		
		}else{
			for(i in 1:num_grps_selected){
				selected_grp_id=which(grps_selected[i]==all_grp_names);
				ctl_name=paste("ModelBuilderTab.group.selector_", selected_grp_id, sep="");
				
				session$userData[["ModelBuilderSets"]][[ctl_name]]=
					sort(union(session$userData[["ModelBuilderSets"]][[ctl_name]], variables_selected));
				updateSelectInput(session, ctl_name, choices=session$userData[["ModelBuilderSets"]][[ctl_name]]);
			}
			
			# Clear group selector
			updateSelectizeInput(session, "ModelBuilderTab.group_selector", selected="");
			
			# Clear variable selector
			updateSelectInput(session, "ModelBuilderTab.available_selector", selected="");
			
		}
	
	});

	#--------------------------------------------------------------------------
	
	#for(i in 1){
	#	grps_ctl_id=paste("ModelBuilderTab.group.name_", i, sep="");
	#	cat("Observing: ", grps_ctl_id, "\n");
	#	observeEvent(input[[grps_ctl_id]], {
	#		cat(grps_ctl_id, ": Group name changed:", input[[grps_ctl_id]], "\n");
	#	});
	#}
	

	observeEvent(input$ModelBuilderTab.group.name_1, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	observeEvent(input$ModelBuilderTab.group.name_2, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	observeEvent(input$ModelBuilderTab.group.name_3, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	observeEvent(input$ModelBuilderTab.group.name_4, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	observeEvent(input$ModelBuilderTab.group.name_5, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	observeEvent(input$ModelBuilderTab.group.name_6, {
		ModelBuilderTab.update_grp_selector(input, session);
	});
	
	
	observeEvent(input$ModelBuilderTab.group.delete_1, {
		ModelBuilderTab.delete_grp_var(input, session, 1);
	});
	observeEvent(input$ModelBuilderTab.group.delete_2, {
		ModelBuilderTab.delete_grp_var(input, session, 2);
	});
	observeEvent(input$ModelBuilderTab.group.delete_3, {
		ModelBuilderTab.delete_grp_var(input, session, 3);
	});
	observeEvent(input$ModelBuilderTab.group.delete_4, {
		ModelBuilderTab.delete_grp_var(input, session, 4);
	});
	observeEvent(input$ModelBuilderTab.group.delete_5, {
		ModelBuilderTab.delete_grp_var(input, session, 5);
	});
	observeEvent(input$ModelBuilderTab.group.delete_6, {
		ModelBuilderTab.delete_grp_var(input, session, 6);
	});
	
	
	observeEvent(input$ModelBuilderTab.group.clear_1, {
		ModelBuilderTab.clear_grp_var(input, session, 1);
	});
	observeEvent(input$ModelBuilderTab.group.clear_2, {
		ModelBuilderTab.clear_grp_var(input, session, 2);
	});
	observeEvent(input$ModelBuilderTab.group.clear_3, {
		ModelBuilderTab.clear_grp_var(input, session, 3);
	});
	observeEvent(input$ModelBuilderTab.group.clear_4, {
		ModelBuilderTab.clear_grp_var(input, session, 4);
	});
	observeEvent(input$ModelBuilderTab.group.clear_5, {
		ModelBuilderTab.clear_grp_var(input, session, 5);
	});
	observeEvent(input$ModelBuilderTab.group.clear_6, {
		ModelBuilderTab.clear_grp_var(input, session, 6);
	});
	
	
	observeEvent(input$ModelBuilderTab.group.copyto_1, {
		ModelBuilderTab.copy_to_dialog(input, session, 1);
		ModelBuilderTab.source_copy_grp_id<<-1;
	});
	observeEvent(input$ModelBuilderTab.group.copyto_2, {
		ModelBuilderTab.copy_to_dialog(input, session, 2);
		ModelBuilderTab.source_copy_grp_id<<-2;
	});
	observeEvent(input$ModelBuilderTab.group.copyto_3, {
		ModelBuilderTab.copy_to_dialog(input, session, 3);
		ModelBuilderTab.source_copy_grp_id<<-3;
	});
	observeEvent(input$ModelBuilderTab.group.copyto_4, {
		ModelBuilderTab.copy_to_dialog(input, session, 4);
		ModelBuilderTab.source_copy_grp_id<<-4;
	});
	observeEvent(input$ModelBuilderTab.group.copyto_5, {
		ModelBuilderTab.copy_to_dialog(input, session, 5);
		ModelBuilderTab.source_copy_grp_id<<-5;
	});
	observeEvent(input$ModelBuilderTab.group.copyto_6, {
		ModelBuilderTab.copy_to_dialog(input, session, 6);
		ModelBuilderTab.source_copy_grp_id<<-6;
	});
	
	
	observeEvent(input$ModelBuilderTab.copy_to_dialog.okButton, {
		removeModal();
		cat("Source Group:", ModelBuilderTab.source_copy_grp_id, "\n");
		cat("Destinations:\n");
		target_groups=input$ModelBuilderTab.copy_to_dialog.group_select;
		print(target_groups);
		
		var_sel_from_src_grp=input[[paste("ModelBuilderTab.group.selector_", ModelBuilderTab.source_copy_grp_id, sep="")]];
		
		group_names=ModelBuilderTab.get_all_group_names(input);
		cat("All group names:\n");
		print(group_names)
		for(grp_id in 1:ModelBuilderTab.num_groups){
			if(any(group_names[grp_id]==target_groups)){
				ctl_name=paste("ModelBuilderTab.group.selector_", grp_id, sep="");
				
				session$userData[["ModelBuilderSets"]][[ctl_name]]=
					sort(union(session$userData[["ModelBuilderSets"]][[ctl_name]], var_sel_from_src_grp));
				updateSelectInput(session, ctl_name, choices=session$userData[["ModelBuilderSets"]][[ctl_name]]);
			}
		}
		
		
	});
	observeEvent(input$ModelBuilderTab.copy_to_dialog.cancelButton, {
		removeModal();
	});
	
	#--------------------------------------------------------------------------
	# Help links
	observeEvent(input$ModelBuilderTab.MainHelp, {
		showModal(
			modalDialog(ModelBuilderTab.MainHelpTxt, title="Model Builder", easyClose=T, size="s")
		);
	});
	observeEvent(input$ModelBuilderTab.ExcludedHelp, {
		showModal(
			modalDialog(ModelBuilderTab.ExcludedHelpTxt, title="Excluded Variables", easyClose=T, size="s")
		);
	});
	observeEvent(input$ModelBuilderTab.AvailableHelp, {
		showModal(
			modalDialog(ModelBuilderTab.AvailableHelpTxt, title="Available Variables", easyClose=T, size="s")
		);
	});
	observeEvent(input$ModelBuilderTab.CovariatesHelp, {
		showModal(
			modalDialog(ModelBuilderTab.CovariatesHelpTxt, title="Covariates", easyClose=T, size="s")
		);
	});
	observeEvent(input$ModelBuilderTab.GroupHelp, {
		showModal(
			modalDialog(ModelBuilderTab.GroupHelpTxt, title="Grouped Variables", easyClose=T, size="m")
		);
	});
}


###############################################################################

ModelBuilderTab.MainHelpTxt=tagList(
		tags$p("Use the Model Builder application to assign your variables into multiple types.",
		"This will help you organize the observations you have collected or the factors in",
		"your designed experiment, for statistical modeling."),
		tags$p("Try to divide your variables into 1 of 4 categories: Excluded, Available, Covariates, or Grouped"),
		tags$p("Click on the help links next to each variable selector control for additional",
		"information on that variable category.")
	);
ModelBuilderTab.ExcludedHelpTxt=tagList(
		tags$p("Excluded variables are columns in your metadata that you do not want to include",
		"in your analyses. Note that variables describing the nature of repeated measures, e.g. offsets,",
		"subject IDs, should be excluded since their utilization is analysis specific. The following",
		"are examples of excludable columns:"),
		tags$ul(
			tags$li("Subject IDs"),
			tags$li("Sample IDs"),
			tags$li("Time Offsets"), 
			tags$li("Dates"),
			tags$li("Raw variables, where the transformed version is preferred"),
			tags$li("Variables that need more curation"),
			tags$li("Variables have too many NAs"),
			tags$li("Variables that you want to save for another analysis")
		),
		tags$p("You can use the arrow buttons to move an excluded variable into or out of the Available",
			"category.  Excluded variables will be removed from all Groups if reclassified as such.")
	);
ModelBuilderTab.AvailableHelpTxt=tagList(
		tags$p("Available variables are open for assignment. All variables will start here."),
		tags$p("Available variables may be assigned Excluded or Covariates, at which time ",
		"they would moved out of availability."),
		tags$p("Since a variable may belong in multiple Grouped sets, assigning a variable into ",
		"a Group, will still leave it available for additonal assigments.")
	);
ModelBuilderTab.CovariatesHelpTxt=tagList(
		tags$p("Covariates are variables that you want to consider immutable with respect ",
		"to your experiment/study but you may need to included in your model to control for. ",
		"For example:"),
		tags$ul(
			tags$li("Sex"),
			tags$li("Age"),
			tags$li("Ethnicity"),
			tags$li("Smoking Status"),
			tags$li("Experimental Treatments"),
			tags$li("Medications"),
			tags$li("Diet")
		),
		tags$p("These are general suggestions, so of course there are always exceptions, depending on",
		"the nature of the hypothesis."),
		tags$p("When a variable is assigned as a Covariate, it is removed from any Group assignment"),
		tags$p("Covariates are always considered predictors, i.e. x's in a linear model")			
	);
ModelBuilderTab.GroupHelpTxt=tagList(
		tags$p("Grouped variables are variables that will be treated as a multivariate response (Y) ",
		"or multiple regression predictors (X).  The goal of assigning variables as Grouped is to be able ",
		"to generate  bidirectional models, where the group variables are fit as predictors in one ",
		"model then fit again as a response in the second model, so their relationship (as a stronger ",
		"predictor or responder) can be ascertained.  In both models, the Covariates are predictors,",
		"since they should be considered immutable by any of the variables in the Group."),
		tags$p("If you have a few available variables, i.e. < ~10, you can move all variables into a single group.",
		"However, if you have a set of variables that appear to be related, then they would benefit from grouping."),
		tags$p("Since grouped variables tend to be correlated, isolating them may help you identify variables",
		"that can be removed without affecting the information represented by that group, through:"),
		tags$ul(
			tags$li("Expert-based variable removal, e.g. knowing one measurement is ", 
				"typically more reliable, or clinically relevant, than another"),
			tags$li("Variable selection through penalized regression"),
			tags$li("Principal Components Analysis (PCA)"),
			tags$li("Missing values tend to be correlated, so it is possible to maximize the sample",
				"size used in analyzing each group separately at first")
		),
		tags$p("Some examples where variables should be grouped togetehr include measurements from:"),
		tags$ul(
			tags$li("A assay panel, e.g. cytokines"),
			tags$li("Related physiology, e.g. pulmonary tests, cardiovascular, blood test"),
			tags$li("A metabolic pathway")
		)
	);

ModelBuilder.available_variables=c();
ModelBuilderTab.num_groups=6;
ModelBuilderTab.source_copy_grp_id=0;

if(!exists("integration")){

	avail_variables=sort(c(
		"Apples",
		"Oranges",
		"SubjectIDs",
		"Pears",
		"SampleType",
		"PrePost",
		"Dates",
		"Times"
	));


	ui = fluidPage(
		mainPanel(
			tabsetPanel(
				tabPanel("Data", ""),
				tabPanel("Curation", ""),
				tabPanel("Study", ""),
				ModelBuilderTab(avail_variables)
			)
		)
	);

	ModelBuilderTab.get_data_colnames=function(){
		return(colnames(test.matrix));
	}

	###############################################################################

	server = function(input, output, session) {

		ModelBuilderTab.sets<<-list();	
		ModelBuilderTab.sets[["ModelBuilderTab.available_selector"]]<<-avail_variables;

		observe_ModelBuilderTabEvents(input, output, session);
	}

	###############################################################################
	# Launch
	shinyApp(ui, server);
}
