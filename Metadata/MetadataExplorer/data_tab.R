library(shiny)
library(DT)

DataTab=function(){

	tabPanel("Data",
		sidebarLayout(
			sidebarPanel(
			
				tags$b("Samples x Variables: "),
				textOutput("DataTab.rowxcol"),
				tags$br(),
				
				radioButtons(
					inputId="DataTab.col_order_radiobutton",
					label="Variable Order:",
					choices=c("Original", "Alphabetic", "Non-NA", "Entropy"),
					selected="Original"
				),
				
				selectInput(
					"DataTab.disp_col_selector", 
					"Variables Shown:",
					choices=colnames(data_matrix),
					multiple=T, selectize=F, size=20,
					selected=colnames(data_matrix)[1:10]),
					
				width=2
				),
			mainPanel(
				tags$br(),
				dataTableOutput("DataTab.table")
			)
		)	
	)
}

observe_DataTabEvents=function(input, output, session){

	output$DataTab.table=renderDT({
		DT::datatable(data_matrix[, input$DataTab.disp_col_selector, drop=F])
	});
	
	output$DataTab.rowxcol=renderText({ paste(nrow(data_matrix), " x ", ncol(data_matrix), sep="")});
	
	observeEvent(input$DataTab.col_order_radiobutton,{
		selected=input$DataTab.disp_col_selector;
		DataTab.column_order<<-DataTab.get_column_order(data_matrix, input$DataTab.col_order_radiobutton);
		updateSelectInput(session, "DataTab.disp_col_selector", choices=DataTab.column_order, selected=selected);
	});
	
}

###############################################################################

DataTab.get_column_order=function(mat, ordering){
	varnames=colnames(mat);
	if(ordering=="Original"){
		var_ordering=varnames;
	}else if(ordering=="Alphabetic"){
		var_ordering=sort(varnames);
	}else if(ordering=="Non-NA"){
		var_ordering=varnames[order(decreasing=F, get_numNAs(mat))];
	}else if(ordering=="Entropy"){
		var_ordering=varnames[order(decreasing=T, get_entropy(mat))];
	}else{
		var_ordering=varnames;
	}

	print(var_ordering);
	return(var_ordering)
}

#------------------------------------------------------------------------------

get_entropy=function(mat){

	ncolumns=ncol(mat);
	entropy=numeric(ncolumns);
	names(entropy)=colnames(mat);
	
	for(i in 1:ncolumns){
		cur_col=mat[,i];
		
		asnumeric=as.numeric(cur_col);
		nas_ix=is.na(asnumeric);
		num_numeric=sum(!nas_ix);
		
		# Binning strategy
		if(num_numeric>0){
			non_nas_vals=asnumeric[!nas_ix];
		}else{
			non_nas_vals=cur_col[!is.na(cur_col)];
			non_nas_vals=as.numeric(as.factor(non_nas_vals));
		}
		
		ranges=range(non_nas_vals);
		bins=seq(ranges[1], ranges[2], length.out=20);
		hist_rec=hist(non_nas_vals, breaks=bins, plot=F);
		
		totals=sum(hist_rec$counts);
		prob=hist_rec$counts/totals;
		prob=prob[prob>0];
		entropy[i]=-sum(log2(prob)*prob);
	}
	
	return(entropy);
}

#------------------------------------------------------------------------------

get_numNAs=function(mat){
	num_nas=apply(mat, 2, function(x){ sum(is.na(x))});
	print(num_nas);
	return(num_nas);
}

###############################################################################

ui = fluidPage(
	mainPanel(
		tabsetPanel(
			DataTab(),
			tabPanel("Curation", ""),
			tabPanel("Study",""),
			tabPanel("Model Explorer", "")
		)
	)
);

test.cols=c("SampleID", "SubjectID", "Time", "Age", "Sex", 
			"PercPred", "Probability", "Variable1", "Systolic", "Diastolic", 
			"Diet", "Occupation");
nts=50;

test.matrix=cbind(
	paste("samp", rpois(nts, 100000), sep=""),
	paste("subj", rpois(nts, 10), sep=""),
	rpois(nts, 30),
	round(rnorm(nts, 40),2),
	rbinom(nts, 1, .6),
	
	round(rnorm(nts, 100), 2),
	round(runif(nts, 0,1),4),
	round(rnorm(nts, 85, 1),2),
	round(rnorm(nts, 120, 20),1),
	round(rnorm(nts, 80, 20),1),
	
	sample(c("Apples", "Bananas", "Cucumbers", NA), nts, replace=T),
	sample(c("Fireman", "Law Enforcement", "Dentist", "Teacher", "Nuclear Scientist", "Lawyer", "Politician", NA), nts, replace=T)
);

colnames(test.matrix)=test.cols;
test.matrix=test.matrix[order(test.matrix[,"SubjectID"]),];

#print(test.matrix);
data_matrix=test.matrix;
DataTab.column_order=colnames(data_matrix);

###############################################################################

server = function(input, output, session) {
	observe_DataTabEvents(input, output, session);
}

###############################################################################
# Launch
shinyApp(ui, server);
