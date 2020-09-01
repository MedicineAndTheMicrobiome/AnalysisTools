rm(list=ls(all=TRUE));

library(shiny);
integration=TRUE;

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\import_export_tab.R");
cat("ImportExportTab Sourced.\n");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\data_tab.R");
cat("DataTab Sourced.\n");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_tab.R");
cat("CurationTab Sourced.\n");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\study_tab.R");
cat("StudyTab Sourced.\n");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\model_builder_tab.R");
cat("ModelBuilderTab Sourced.\n");

###############################################################################

ui = fluidPage(
	mainPanel(
		tabsetPanel(
			ImportExportTab(),
			DataTab(),
			CurationTab(),
			StudyTab(),
			ModelBuilderTab()
		)
	)
);

###############################################################################

server = function(input, output, session) {

	session$userData[["ModelBuilderSets"]]=ModelBuilderTab.init_variables();
	session$userData[["Metadata"]]=matrix(NA, ncol=0, nrow=0);

	observe_ImportExportTabEvents(input, output, session);
	observe_DataTabEvents(input, output, session);
	observe_CurationTabEvents(input, output, session);
	observe_StudyTabEvents(input, output, session);
	observe_ModelBuilderTabEvents(input, output, session);
}

###############################################################################
# Launch
shinyApp(ui, server);
