library(shiny)

fluidPage(

  titlePanel("Internal Tools"),

  sidebarPanel(

	fileInput("file1", "Choose 1st TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file2", "Choose 2nd TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file3", "Choose 3rd TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file4", "Choose 4th TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file5", "Choose 5th TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file6", "Choose 6th TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file7", "Choose 7th TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv")),
	fileInput("file8", "Choose 8th TSV File",
                multiple = FALSE,
                accept = c("text/tsv",
                         ".tsv"))
  ),

  mainPanel(
	tags$p("When the tool has completed running, you may download the output:"),
	downloadButton("downloadButton", "Download")
  )
)