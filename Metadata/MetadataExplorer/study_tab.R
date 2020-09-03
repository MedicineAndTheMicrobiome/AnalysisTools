
library(shiny)


StudyTab=function(){

	tabPanel("Study",
		sidebarLayout(
			sidebarPanel(
				tags$h1(),
				radioButtons("StudyTab.StudyType_radioButton",
					label="Please identify your study type:",
					choices=c("Cross-Sectional", "Longitudinal", "Paired"),
					selected=NULL),
				uiOutput("StudyTab.help_text_ui")
			),
			mainPanel(
				uiOutput("StudyTab.study_specific_ui")
			)
		)	
	)

}

observe_StudyTabEvents=function(input, output, session){

	observeEvent(input$StudyTab.StudyType_radioButton,{
		cat("StudyType Radio Button Clicked.\n");
		
		# Update the help text
		output$StudyTab.help_text_ui=renderUI(
			switch(input$StudyTab.StudyType_radioButton,
				"Cross-Sectional" = 
						{tags$html(StudyTab.crosssection_help_text);},
				"Longitudinal" =
						{tags$html(StudyTab.longitudinal_help_text);},
				"Paired" = 
						{tags$html(StudyTab.paired_help_text);}
			)
		);
		
		# Reset the prior values
		switch(input$StudyTab.StudyType_radioButton,
				"Cross-Sectional" = {
					updateSelectInput(session, "StudyTab.StudyType_subjectID.cs", 
							selected=input$StudyTab.StudyType_subjectID.cs);},
				"Longitudinal" = {
					updateSelectInput(session, "StudyTab.StudyType_subjectID.long", 
							selected=input$StudyTab.StudyType_subjectID.long);
					updateSelectInput(session, "StudyTab.StudyType_TimeOffsets", 
							selected=input$StudyTab.StudyType_TimeOffsets);
				},
				"Paired" = {
					updateSelectInput(session, "StudyTab.StudyType_subjectID.paired", 
							selected=input$StudyTab.StudyType_subjectID.paired);
					updateSelectInput(session, "StudyTab.StudyType_Pairings", 
							selected=input$StudyTab.StudyType_Pairings);
				}
			);		
		
		# (Re)draw Select variable controls
		output$StudyTab.study_specific_ui=renderUI({
			
			avail_variables=StudyTab.get_data_colnames(session);
			
			fluidPage(
			tags$h1(),
			tags$h4("Identify the necessary variables from your dataset:"),
		
			switch(input$StudyTab.StudyType_radioButton,
				"Cross-Sectional" = 
						{
						selectInput("StudyTab.StudyType_subjectID.cs", "Subject IDs", choices=avail_variables)
						},
				"Longitudinal" =
						{
						tagList(
						selectInput("StudyTab.StudyType_subjectID.long", "Subject IDs", choices=avail_variables),
						selectInput(
							"StudyTab.StudyType_TimeOffsets", "Time Offets", choices=avail_variables)
						)
						},
				"Paired" = 
						{
						tagList(
						selectInput("StudyTab.StudyType_subjectID.paired", "Subject IDs", choices=avail_variables),
						selectInput(
							"StudyTab.StudyType_Pairings", "Pairing Criteria", choices=avail_variables)
						)
						}
				),
				
			tags$h1(),
			tags$h5("Let us check your selected columns for suitability."),
			actionButton("StudyTab.StudyType_checkButton", "Check"),
			tags$h1(),
			wellPanel(
				uiOutput("StudyTab.StudyType_checkResults_ui")
			)
							
			);
			
		})
		
		# Clear the output, just in case different radio button was selected
		output$StudyTab.StudyType_checkResults_ui=renderUI(tagList());
		
		
	});

	observeEvent(input$StudyTab.StudyType_checkButton,{
	
		session$userData[["StudyType"]][["Type"]]=input$StudyTab.StudyType_radioButton;
			
		switch(input$StudyTab.StudyType_radioButton,
				"Cross-Sectional" = 
						{
						
						sbj_ids=StudyTab.get_column(input$StudyTab.StudyType_subjectID.cs, session);
										
						output$StudyTab.StudyType_checkResults_ui=
							renderUI(tagList(
							
									tags$h3("Subject IDs"),
									fluidRow(
									column(3,
										selectInput("StudyTab.exampleSubjectIDs_select", 
											"Examples: ", 
											choices=head(sbj_ids,1000), selectize=F, size=15)),
									column(9, plotOutput("StudyTab.subjectid_validate_plots"))
									)));
									
						
		
						output$StudyTab.subjectid_validate_plots=renderPlot({
								StudyTab.subjectid_validate(sbj_ids, "cross-sectional");
							});
													
						session$userData[["StudyType"]][["SampleID.cs"]]=input$StudyTab.StudyType_subjectID.cs;
						},
				"Longitudinal" =
						{
						sbj_ids=StudyTab.get_column(input$StudyTab.StudyType_subjectID.long, session);
						time_offsets=StudyTab.get_column(input$StudyTab.StudyType_TimeOffsets, session);
						
						output$StudyTab.StudyType_checkResults_ui=
							renderUI(tagList(
							
									tags$h3("Subject IDs"),		
									fluidRow(
									column(3,
										selectInput("StudyTab.exampleSubjectIDs_select", 
											"Examples: ", 
											choices=head(sbj_ids,1000), selectize=F, size=17)),
									column(9,
										plotOutput("StudyTab.subjectid_validate_plots")),
									),
									
									tags$h3("Time Offsets"),
									fluidRow(
									column(3,
										selectInput("StudyTab.exampleTimeOffsets_select", 
											"Examples: ", 
											choices=head(time_offsets,1000), selectize=F, size=17)),
											
									column(9, plotOutput("StudyTab.timeoffset_validate_plots"))
									)
									
									));
									
						output$StudyTab.subjectid_validate_plots=renderPlot({
							StudyTab.subjectid_validate(sbj_ids, "longitudinal");
							});						
						
						output$StudyTab.timeoffset_validate_plots=renderPlot({
								StudyTab.timeoffset_validate(sbj_ids, time_offsets);
							});
									
						session$userData[["StudyType"]][["SampleID.long"]]=input$StudyTab.StudyType_subjectID.long;
						session$userData[["StudyType"]][["TimeOffsets"]]=input$StudyTab.StudyType_radioButton;	
											
						},
				"Paired" = 
						{
						sbj_ids=StudyTab.get_column(input$StudyTab.StudyType_subjectID.paired, session);
						pairings=StudyTab.get_column(input$StudyTab.StudyType_Pairings, session);
						
						output$StudyTab.StudyType_checkResults_ui=
							renderUI(tagList(
							
									tags$h3("Subject IDs"),		
									fluidRow(
									column(3,
										selectInput("StudyTab.exampleSubjectIDs_select", 
											"Examples: ", 
											choices=head(sbj_ids,1000), selectize=F, size=17)),
									column(9,
										plotOutput("StudyTab.subjectid_validate_plots")),
									),
									
									tags$h3("Pairings"),
									fluidRow(
									column(3,
										selectInput("StudyTab.examplePairings_select", 
											"Examples: ", 
											choices=head(pairings,1000), selectize=F, size=17)),
											
									column(9, plotOutput("StudyTab.pairings_validate_plots"))
									)
									
									));
									
						output$StudyTab.subjectid_validate_plots=renderPlot({
							StudyTab.subjectid_validate(sbj_ids, "paired");
							});						
						
						output$StudyTab.pairings_validate_plots=renderPlot({
								StudyTab.pairings_validate(sbj_ids, pairings);
							});
							
						session$userData[["StudyType"]][["PairingCriteria"]]=input$StudyTab.StudyType_radioButton;
						session$userData[["StudyType"]][["SampleID.paired"]]=input$StudyTab.StudyType_subjectID.paired;

						}
				)
	});
}


###############################################################################

StudyTab.crosssection_help_text=
	tags$html(tags$b("Cross-sectional study:"), tags$p("A study assuming a single observation per subject. "), 
		tags$p("This type of study must have a unique subject ID for each sample(row)."));
StudyTab.longitudinal_help_text=
	tags$html(tags$b("Longitudinal study:"), tags$p("A study with multiple repeated measures for each subject,",
		" collected over time."), tags$p("This type of study will have multiple subject IDs, but each subject will ",
		" have multiple times (offsets) associated with each sample(row)."));
StudyTab.paired_help_text=
	tags$html(tags$b("Paired study:"), tags$p("A study with 2 measures for each subject. ",
		"One measurement type is compared against the other from the same individual."), 
		tags$p("This type of study will have 2 samples(rows) per subject, and a pairing criteria ",
		" that differentiates the 2 samples from each other."));



StudyTab.get_column=function(var_name, session){
	return(session$userData[["Metadata"]][,var_name]);
}	
StudyTab.get_data_colnames=function(session){
	return(colnames(session$userData[["Metadata"]]));
}




#------------------------------------------------------------------------------

StudyTab.subjectid_validate=function(subject_ids, study_type){
	
	na_ix=is.na(subject_ids);
	num_nas=sum(na_ix);
	total_samples=length(subject_ids);
	samp_nonas=subject_ids[!na_ix];
	num_nonas=length(samp_nonas);
	uniq_samp=unique(samp_nonas);
	num_uniq=length(uniq_samp);
	num_dups=num_nonas-num_uniq;
	
	bar_values=c(total_samples, num_nonas, num_uniq, num_dups, num_nas);
	max_bar_value=max(bar_values);
	prop_uniq=num_uniq/num_nonas;
	
	par(mar=c(7,4,4,1));
	
	if(study_type=="cross-sectional"){
		
		mids=barplot(
			bar_values,
			names.arg=c("Rows", "Non-NAs", "Unique", "Duplicates", "NAs"),
			col=c("grey", "grey", "blue", "red", "grey"),
			ylab="Counts",
			ylim=c(0, max_bar_value*1.1)
		);
		
		text(mids, bar_values, bar_values, pos=3);
		
		mtext("In a cross-sectional study, all of your non-NA Subject IDs should be unique.",
			side=1,
			line=3,
			);
		
		text_col=ifelse(prop_uniq<0.90, "red", "blue");
		mtext(paste(round(100*prop_uniq,2), "% of your Subject IDs are unique.", sep=""),
			side=1,
			line=5,
			col=text_col
		);
	
	}else if(study_type=="longitudinal"){
		
		mids=barplot(
			bar_values,
			names.arg=c("Rows", "Non-NAs", "Unique", "Duplicates", "NAs"),
			col=c("grey", "grey", "red", "blue", "grey"),
			ylab="Counts",
			ylim=c(0, max_bar_value*1.1)
		);
		
		text(mids, bar_values, bar_values, pos=3);
		
		mtext("In a longitudinal study, you would expect to have more than 1 measurement per subject.",
			side=1,
			line=3,
			);
			
		text_col=ifelse(prop_uniq>0.10, "red", "blue");
		mtext(paste(round(100*prop_uniq,2), "% of your Subject IDs do not have repeated measurements.", sep=""),
			side=1,
			line=5,
			col=text_col
		);
		
	}else if(study_type=="paired"){
	
		mids=barplot(
			bar_values,
			names.arg=c("Rows", "Non-NAs", "Unique", "Duplicates", "NAs"),
			col=c("grey", "grey", "red", "blue", "grey"),
			ylab="Counts",
			ylim=c(0, max_bar_value*1.1)
		);
		
		text(mids, bar_values, bar_values, pos=3);
		
		mtext("In a paired study, most of your non-NA Subject IDs should have a matching pair.",
			side=1,
			line=3,
			);
			
		text_col=ifelse(prop_uniq>0.10, "red", "blue");
		mtext(paste(round(100*prop_uniq,2), "% of your Subject IDs are not paired.", sep=""),
			side=1,
			line=5,
			col=text_col
		);
	}

}

StudyTab.timeoffset_validate=function(subject_ids, timeoffsets){
	
	is_numeric=is.numeric(timeoffsets);
	offsets_val=as.numeric(timeoffsets);
	na_ix=is.na(offsets_val);
	nonna_offsets_val=offsets_val[!na_ix];
	
	num_nonna=length(nonna_offsets_val);
	
	if(num_nonna>0 && is_numeric){
		par(mfrow=c(2,1));
		
		par(mar=c(4,4,2,1));
		hist(offsets_val, xlab="Time Offsets", main="", ylab="Counts");
		
		par(mar=c(8,4,2,1));
		hist(table(subject_ids), xlab="Measures Per Subject", ylab="Counts", main="")
		
		mtext("The distribution of time offsets and measurements per subject", side=1, line=5, col="blue");
		mtext("will depend on your sampling frequency and rate of subject drop out.", side=1, line=6, col="blue");
		
		
	}else{
		cat_counts=table(timeoffsets);
		plot(0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="", bty="n",
			xlim=c(-1,1), ylim=c(-1,1));
		text(0,0, "Error: Time Offsets are not numeric.", font=2, col="red");
	}
}

StudyTab.pairings_validate=function(sbj_ids, pairings){
	
	nonna_sbj_ids=!is.na(sbj_ids);
	uniq_subj_ids=unique(sbj_ids[nonna_sbj_ids]);
	
	nonna_pairings=!is.na(pairings);
	uniq_pairings=unique(pairings[nonna_pairings]);
	num_uniq_pairings=length(uniq_pairings);
	
	print(uniq_pairings);
	
	if(num_uniq_pairings!=2){
		plot(0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="", bty="n",
			xlim=c(-1,1), ylim=c(-1,1));
		text(0,0, "Error: There are not exactly 2 pairing types.", font=2, col="red");
	}else{
		pairing_list_cts=list();
		pairing_list_cts[[uniq_pairings[1]]]=0;
		pairing_list_cts[[uniq_pairings[2]]]=0;
		pairing_list_cts[["Both"]]=0;
		
		error=F;
		for(sbj_ix in uniq_subj_ids){
		
			ix=(sbj_ix == sbj_ids);
			matched=pairings[ix];
			matched=matched[!is.na(matched)];
			num_types=length(matched);
			
			if(num_types==1){
				pairing_list_cts[[matched]]=pairing_list_cts[[matched]]+1;
			}else if(num_types==2){
				if(matched[1]!=matched[2]){
					pairing_list_cts[["Both"]]=pairing_list_cts[["Both"]]+1;
				}else{
					error=T;
					break;
				}
			}else if(num_types==0){
			
			}else if(num_types>2){
				error=T;
				break;
			}
		}
		
		
		if(error){
			plot(0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="", bty="n",
				xlim=c(-1,1), ylim=c(-1,1));
			text(0,0, "Error: Sample(s) with more than one potential pairing.", font=2, col="red");
		}else{
		
			bar_values=c(
				pairing_list_cts[[uniq_pairings[1]]],
				pairing_list_cts[["Both"]],
				pairing_list_cts[[uniq_pairings[2]]]
				);
				
			par(mar=c(7,4,1,1));
			mids=barplot(
				bar_values,
				names.arg=c(
				paste("Only:\n", uniq_pairings[1], sep=""), 
				"Complete\nPairs",
				paste("Only:\n", uniq_pairings[2], sep="")),
				col=c("red", "purple", "blue"),
				ylim=c(0, max(bar_values)*1.1)
			);
			
			text(mids, bar_values, bar_values, pos=3);
			
			prop_paired=pairing_list_cts[["Both"]]/(pairing_list_cts[[uniq_pairings[1]]]+pairing_list_cts[[uniq_pairings[2]]]);
			msg_col=ifelse(prop_paired<.9,	"red", "blue");
				
			mtext(paste(round(100*prop_paired,2), "% of subjects are paired.", sep=""), side=1, line=3, col=msg_col);

		}
	}
	
}

	
if(!exists("integration")){
	ui = fluidPage(
		mainPanel(
			tabsetPanel(
				tabPanel("Data", ""),
				tabPanel("Curation", ""),
				StudyTab(),
				tabPanel("Model Explorer", "")
			)
		)
	);

	avail_variables=c(
		"Apples",
		"Oranges",
		"SubjectIDs",
		"Pears",
		"SampleType",
		"PrePost",
		"Dates",
		"Times"
	);
	
	test.nrows=200;

	test.cs_sbjIds={a=.05; sample(c(rep(NA,test.nrows*a), paste("s_", 1:(test.nrows*(1-a)), sep="")), test.nrows)};
	test.long_sbjIds=sample(test.cs_sbjIds[1:20], test.nrows, replace=T);
	test.offsets=sample({a=.7;as.integer(c(rexp(test.nrows*a,.0005)*.005,rpois(test.nrows*(1-a),35)))}, test.nrows);

	test.pairings_sbjIds=paste("ps_", rep(1:(test.nrows/2),2), sep="");
	test.pairings_criteria={a=.20; b=.15; 
		c(
			rep(NA, test.nrows/2*a),
			rep("Before", test.nrows/2*(1-a)),
			rep("After", test.nrows/2*(1-b)),
			rep(NA, test.nrows/2*b)
		)};
		
	test.matrix=cbind(test.cs_sbjIds, test.pairings_sbjIds, test.pairings_criteria, test.long_sbjIds, test.offsets);
	print(test.matrix);

	avail_variables=c(
		StudyTab.get_data_colnames(),
		"Apples",
		"Oranges",
		"SubjectIDs",
		"Pears",
		"SampleType",
		"PrePost",
		"Dates",
		"Times"
	);

	###############################################################################

	server = function(input, output, session) {

		observe_StudyTabEvents(input, output, session, avail_variables);
	}

	###############################################################################
	# Launch
	shinyApp(ui, server);

}