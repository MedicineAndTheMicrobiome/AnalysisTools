library(shiny)

fluidPage(

  titlePanel("Internal Tools"),

  sidebarPanel(

	tags$h5(tags$b("Select your input files:")),
	tags$h6(tags$p("Use the checkbox to disable an uploaded file.")),
	tags$h6(tags$p("You may replaced a previously uploaded file by selecting a new file to replace it.")), 

	fluidRow(
		column(1, checkboxInput("file1use_cb", "", value=F)),
		column(8, fileInput("file1_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file2use_cb", "", value=F)),
		column(8, fileInput("file2_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file3use_cb", "", value=F)),
		column(8, fileInput("file3_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file4use_cb", "", value=F)),
		column(8, fileInput("file4_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file5use_cb", "", value=F)),
		column(8, fileInput("file5_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file6use_cb", "", value=F)),
		column(8, fileInput("file6_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file7use_cb", "", value=F)),
		column(8, fileInput("file7_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),

	fluidRow(
		column(1, checkboxInput("file8use_cb", "", value=F)),
		column(8, fileInput("file8_upld", NULL, multiple = FALSE, accept = c("text/tsv", ".tsv")))
		),
  ),

  mainPanel(
	tags$p("When you are ready to run, press this button:"),
	fluidRow(
		column(2, actionButton("runButton", "Run")),
		column(1, span(textOutput("runNotReadyTxt"), style="color:red")),
		column(1, span(textOutput("runReadyTxt"), style="color:green"))
	),
	tags$p(HTML("&nbsp")),
	
	tags$p("When the tool has completed running, you may download the output:"),
	fluidRow(
		column(2, downloadButton("downloadButton", "Download")),
		column(1, span(textOutput("downloadNotReadyTxt"), style="color:red")),
		column(1, span(textOutput("downloadReadyTxt"), style="color:green"))
	),
	tags$p(HTML("&nbsp")),
	
	tags$p("Standard Output from Running Tool:"),
	verbatimTextOutput("appstdoutReadyTxt", placeholder=T)
  )
)