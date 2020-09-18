
library(shiny)

###############################################################################

DataTransformationDialogBoxUI=function(id, in_values, in_varname, default_trans="x"){
	
	ns=NS(id);
	
	modalDialog(
		tagList(
			 titlePanel(paste("Data Transformation: ", in_varname, sep="")),
			 helpText("Please find a suitable transformation for your variable:"),
			 fluidRow(
				column(5,
					wellPanel(
						renderPlot({hist(in_values, main="Original", xlab=in_varname);}),
						htmlOutput(ns("DTDB.untransformed_shapirowilks_pval"))
					)
				),
				column(2,
					verticalLayout(
						tags$h1(""),tags$h1(""),tags$h1(""),
						selectInput(ns("DTDB.transformation_select"), 
								  label = "Transformations:",
								  choices = DTDB.transformations,
								  selected = default_trans,
								  size=8, selectize=F
								  )
					)
				),
				column(5,
					wellPanel(
						plotOutput(ns("DTDB.transformed_histogram")),
						htmlOutput(ns("DTDB.transformed_shapirowilks_pval"))
					)
				)
			)
		),
		footer=fluidRow(
				column(1, actionButton(ns("DTDB.helpButton"), label="Help")),
				column(7),
				column(2, actionButton(ns("DTDB.saveAsButton"), label="Save As...")),
				column(2, actionButton(ns("DTDB.cancelButton"), label="Cancel"))
			),
		size="l"
	);
	
}

###############################################################################

DataTransformationDialogBoxServer=function(id, invalname, invarname, transformname, save_transf_data_callback){

	moduleServer(
	
		id,
		
		function(input, output, session){
		
			ns=NS(id);
		
			observeEvent(input$DTDB.helpButton, {
				cat("Dialog Help Pushed.\n");
				removeModal();
				showModal(modalDialog(title="Data Transformation", 
					DTDB.help_txt, 
					footer=actionButton(ns("DTDB.dismissHelp"), label="OK"),
					size="l"
					));
			});
			
			observeEvent(input$DTDB.cancelButton, {
				cat("Dialog Cancel Pushed.\n");
				removeModal();
			});
			
			observeEvent(input$DTDB.saveAsButton, {
				cat("Dialog SaveAs Pushed.\n");
				
				if(!is.null(input$DTDB.transformation_select)){
					removeModal();
					sug_name=DTDB.suggest_name(input$DTDB.transformation_select, session$userData[[invarname]]);
					cat("Suggested Name: ", sug_name, "\n");
					showModal(RenameDialogBoxUI(
						id="data_rename",
						old_name=session$userData[[invarname]],
						suggested_name=sug_name,
						mesg="Please specify a new name for your transformed variable."));
				}else{
					showModal(CantSave_DialogBox("No transformation specified"));
				}
			});
			
			observeEvent(input$DTDB.transformation_select, {
				cat("Transformation selector touched.\n");
				cat(input$DTDB.transformation_select, "\n");
			
				warnings_table=session$userData[["Curation"]][["WarningsTable"]];
			
				session$userData[[transformname]]=input$DTDB.transformation_select;
			
				transformed_values=DTDB.transform(input$DTDB.transformation_select, session$userData[[invalname]]);
				nas=is.na(transformed_values);
				num_nas=sum(nas);
				
				orig_pval=tryCatch({
					res=shapiro.test(session$userData[[invalname]]);
					res$p.value;
				}, warning=function(w){}, error=function(e){});
				
				trans_pval=tryCatch({
					res=shapiro.test(transformed_values);
					res$p.value;
				}, warning=function(w){}, error=function(e){});
				
				bad_trans=F;
				
				if(is.null(trans_pval)){
					trans_pval=0;
					bad_trans=T;
				}
				
				if(is.null(orig_pval)){
					orig_pval=0;
				}
				
				if(orig_pval>.10){
					orig_col="green";
				}else{
					orig_col="red";
				}
				
				if(trans_pval>.10){
					trans_col="green";
				}else{
					trans_col="red";
				}
			
				if(trans_pval>orig_pval){
					trans_bold=c("<b>", "</b>");
					orig_bold=c("","");
				}else{
					orig_bold=c("<b>", "</b>");
					trans_bold=c("","");
				}
			
				cat("invarname: ", session$userData[[invarname]], "\n");
				xlabel=gsub("x", session$userData[[invarname]], input$DTDB.transformation_select);
				
				if(bad_trans){
					output$DTDB.transformed_histogram=renderPlot({
						plot(0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main="", xlim=c(-1,1), ylim=c(-1,1));
						text(0,0, "Transformation not applicable\nfor values in this range.");
						}
					)
				}else{
					output$DTDB.transformed_histogram=renderPlot(hist(transformed_values, xlab=xlabel, main="Transformed"));
				}
				
				output$DTDB.untransformed_shapirowilks_pval=renderText({
					paste(orig_bold[1], "<p style=\"color:", orig_col, "\">Shapiro-Wilks: p-value = ", signif(orig_pval,3),
						"</p>", orig_bold[2], sep="")});
				
				output$DTDB.transformed_shapirowilks_pval=renderText({
					paste(trans_bold[1], "<p style=\"color:", trans_col, "\">Shapiro-Wilks: p-value = ", signif(trans_pval,3),
						"</p>", trans_bold[2], sep="")});
				
			});

			#--------------------------------------------------------------------------
			# Help Events
			
			observeEvent(input$DTDB.dismissHelp, {
				cat("Dismissed Help.\n");
				removeModal();
				showModal(DataTransformationDialogBoxUI(id, session$userData[[invalname]], session$userData[[invarname]]), session);
				updateSelectInput(session, ns("DTDB.transformation_select"), selected=input$DTDB.transformation_select);
			});

			#--------------------------------------------------------------------------
			# CS, Can't save Events
			
			observeEvent(input$CS.okButton, {
					showModal(DataTransformationDialogBoxUI(id, session$userData[[invalname]], session$userData[[invarname]]), session);
			});
		
			
			}
		 );
}

###############################################################################

DTDB.help_txt=tags$html(
		tags$h4("The Goal of Transforming Data"),
		tags$p(
			"Data transformations are crucial for statistical analyses ",
			"because many models depend on normally distributed measurement ",
			"errors, or the relationship between two ",
			"variables, such as a response and a set of predictors, to ",
			"have a linear relationship."
		),
		tags$br(),
		tags$h4("Available Transformations"),
		tags$p(
			tags$b("sqrt(x)"), ": Counts tend to have poisson distributions and may benefit ",
			"from this transformation by variance stabilization.  ie. ",
			"constant variance independent of the value of x.  Measurements of ",
			"area (e.g., L x W, or r^2) also benefit from this transformation."),
		tags$p(
			tags$b("ln(x)"), " and ", tags$b("log10(x)"), ": Useful when a value has an exponential response to an effect, ",
			"or x is a ratio of a/b.  The choice of natural or a base-10 log is not ",
			"crucial, but the transformed values may be more intuitive under some cirumstances ",
			"when base-10 is applied."),
		tags$p(
			tags$b("ln(x+1)"), " and ", tags$b("log10(x+1)"), ": Measurements may have a value of 0",
			"as a result of assay sensitivity, or sampling size.  (Sample size of the",
			"items being measured, not the number of subjects.)  When 0's are present",
			"the log transformations will return -Inf, so adding 1 to all values will",
			"address this issue."),
		tags$p(
			tags$b("logit(x)"), ": Probabilities and proportions may benefit from this transformation,",
			"especially when the mean is not centered around .5.  Range of transformation",
			"can span from -Inf to +Inf."),
		tags$p(
			tags$b("arcsin(sqrt(x))"), ": Similar to logit.  Use at your own discretion. Range of ",
			"transformation can span from 0 to pi."
		),
		tags$br(),
		tags$h4("How to Explore"),
		tags$p(
			"The left and right panels contain a histogram/distribution of the variable's original and ",
			"transformed values, respectively. The Shapiro-Wilks test, a statistical test for normality ",
			"(where if the p-value > 0.05, the distribution cannot be rejected as non normal), is displayed",
			"in each panel.  If the transformation has a greater p-value (more normal), than the ",
			"variable's original value, then the text will be displayed in bold.  If either original or transformed ",
			"variable may be considered normally distributed, each will be colored green, or else red."),
		tags$p(
			"Remember Occam's Razor when chosing to accept a transformation. ",
			"If the original value is more or less normally distributed, then leave it alone.",
			"If the transformed value does not significantly improve it's normality, then accept the original values."
		)
);


DTDB.transformations=c(
	"sqrt(x)",
	"ln(x)",
	"ln(x+1)",
	"log10(x)",
	"log10(x+1)",
	"logit(x)",
	"arcsin(sqrt(x))",
	"x"
);

DTDB.transform=function(trans, x){

	ret=tryCatch({
		if(trans=="sqrt(x)"){
			trans=(sqrt(x));
		}else if(trans=="ln(x)"){
			trans=(log(x));
		}else if(trans=="ln(x+1)"){
			trans=(log(x+1));
		}else if(trans=="log10(x)"){
			trans=(log10(x));
		}else if(trans=="log10(x+1)"){
			trans=(log10(x+1));
		}else if(trans=="logit(x)"){
			trans=(log(x/(1-x)));
		}else if(trans=="arcsin(sqrt(x))"){
			trans=(asin(sqrt(x)));
		}else if(trans=="x"){
			trans=x;
		}
	});

	return(ret);

}

DTDB.suggest_name=function(trans, varname){

	if(trans=="sqrt(x)"){
		sug_name=paste("sqrt_", varname, sep="");
	}else if(trans=="ln(x)"){		
		sug_name=paste("ln_", varname, sep="");
	}else if(trans=="ln(x+1)"){
		sug_name=paste("ln_", varname, "_p1", sep="");
	}else if(trans=="log10(x)"){
		sug_name=paste("log10_", varname, sep="");
	}else if(trans=="log10(x+1)"){
		sug_name=paste("log10_", varname, "_p1", sep="");
	}else if(trans=="logit(x)"){
		sug_name=paste("logit_", varname, sep="");
	}else if(trans=="arcsin(sqrt(x))"){
		sug_name=paste("asin_sqrt_", varname, sep="");
	}else if(trans=="x"){
		sug_name=varname;
	}

	return(sug_name);

}

if(!exists("integration")){

	###############################################################################
	# Use this function so we can see what is being returned.
	# As we share code, we will look for this call to quickly identify what is being
	# returned by this set of functions.  This code below will not be executed
	# in the main application.

	SetReturn=function(varname, value, session){
		updateTextInput(session, varname, value=value); 
	}

	# This will be provided, for now hard coded for testing
	type="prob";

	if(type=="lognorm"){
		in_values=exp(rnorm(20, 10, 5));
		in_varname="need_log";
		in_default_trans="log(x)";
	}else if(type=="pois"){
		in_values=rpois(100, 10);
		in_varname="need_sqrt";
		in_default_trans="sqrt(x)";
	}else if(type=="prob"){
		in_values=rbeta(300, .35,.25);
		in_varname="beta";
		in_default_trans="logit(x)";
	}else{
		in_values=norm(100, 6, 3);
		in_varname="already_norm";
		in_default_trans="x";
	}

	print(type);
	print(in_varname);


	in_current_variable_names=c("apple_pie", "peach_cobbler", "portuguese_tart", "shoofly_pie");


	###############################################################################
	# Boiling plate for unit testing

	ui = fluidPage(
		actionButton("startButton", label="Start"),
		textInput("Selected_Transformation", label="Selected Transformation:"),
		textInput("Selected_Variable_Name", label="Selected Variable Name:")
	);

	###############################################################################



	server = function(input, output, session) {

		observeEvent(input$startButton,{
			showModal(DataTransformationDialogBox(in_values, in_varname, in_default_trans), session);
		});

		observe_DataTransformationDialogBoxEvents(input, output, session, in_values, in_varname, in_current_variable_names);
	}

	###############################################################################
	# Launch
	shinyApp(ui, server);
}