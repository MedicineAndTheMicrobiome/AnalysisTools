library(shiny)

function(input, output, session) {
	
	output$downloadButton=downloadHandler(
		filename=function(){
			paste("Results.zip");
		},
		content=function(filename){
			
			cat("Write to this file:\n");
			print(filename);
			
			# This is where the uploaded file is on the server
			filelist=c(
				input$file1$datapath,
				input$file2$datapath,
				input$file3$datapath,
				input$file4$datapath,
				input$file5$datapath,
				input$file6$datapath,
				input$file7$datapath,
				input$file8$datapath);
				
			# This is the name of the file the user uploaded (w/o the full path)
			filenames=c(
				input$file1$name,
				input$file2$name,
				input$file3$name,
				input$file4$name,
				input$file5$name,
				input$file6$name,
				input$file7$name,
				input$file8$name
			)
			
			cat("User File Names:\n");
			print(filenames);
			
			cat("Uploaded File Paths:\n");
			print(filelist);
			
			# Write all the files user uploaded, into a table:  1st col "path on server", 2nd col "user file name"
			file_of_paths=tempfile();
			write.table(cbind(filelist, filenames), file_of_paths, sep="\t", row.names=F, col.names=F, quote=F);
			
			# Get path of where we are allowed to write a temp file
			output_root=tempfile();
			
			# Generate system command and then run it.
			mycmd=paste("Rscript www/AccumulateQCAcrossRuns.r -i ", file_of_paths, " -o ", output_root, sep="");
			mydata=system(mycmd, intern=T);
			
			# Debugging output from STDOUT
			print(mydata);
			
			# Generate a list of files 
			first=paste(output_root, ".accumulated.detail.tsv", sep="");
			second=paste(output_root, ".accumulated.simplify.tsv", sep="");
			
			# Zip file together and place results into the filename specified by shiny
			cat("File name to write our output:\n");
			print(c(first,second));
			zip(zipfile=filename, files=c(first, second));
			cat("Done zipping.\n");
		}
	);

}