library(shiny);
library(shinyjs);

setRunReady=function(input, output, session){
	
	select_arr=files_selected(input, output, session);
	
	if(any(select_arr)){
		output$runNotReadyTxt=renderText({""});
		output$runReadyTxt=renderText({"Ready"});
	}else{
		output$runNotReadyTxt=renderText({"Not Ready"});
		output$runReadyTxt=renderText({""});	
	}
}

observeCheckboxEvents=function(input, output, session){
	observeEvent(input$file1use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file1use_cb", value=F);
		}
		setRunReady(input, output, session);
	});
	observeEvent(input$file2use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file2use_cb", value=F);
		}
		setRunReady(input, output, session);
	});	
	observeEvent(input$file3use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file3use_cb", value=F);
		}
		setRunReady(input, output, session);
	});	
	observeEvent(input$file4use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file4use_cb", value=F);
		}
		setRunReady(input, output, session);
	});	
	observeEvent(input$file5use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file5use_cb", value=F);
		}
		setRunReady(input, output, session);
	});	
	observeEvent(input$file6use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file6use_cb", value=F);
		}
		setRunReady(input, output, session);
	});	
	observeEvent(input$file7use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file7use_cb", value=F);
		}
		setRunReady(input, output, session);
	});
	observeEvent(input$file8use_cb, {
		if(ntb(input$file1_upld$datapath)==""){
			updateCheckboxInput(session, "file8use_cb", value=F);
		}
		setRunReady(input, output, session);
	});
}

observeFileUploadsEvents=function(input, output, session){

	observeEvent(input$file1_upld, {
		updateCheckboxInput(session, "file1use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file2_upld, {
		updateCheckboxInput(session, "file2use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file3_upld, {
		updateCheckboxInput(session, "file3use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file4_upld, {
		updateCheckboxInput(session, "file4use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file5_upld, {
		updateCheckboxInput(session, "file5use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file6_upld, {
		updateCheckboxInput(session, "file6use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file7_upld, {
		updateCheckboxInput(session, "file7use_cb", value=T);
		setRunReady(input, output, session);
		});

	observeEvent(input$file8_upld, {
		updateCheckboxInput(session, "file8use_cb", value=T);
		setRunReady(input, output, session);
		});
		
}

ntb=function(x){
	if(is.null(x)){
		return("");
	}else{
		return(x);
	}
}

files_selected=function(input, output, session){

	file_cb=c(
		input$file1use_cb,
		input$file2use_cb,
		input$file3use_cb,
		input$file4use_cb,
		input$file5use_cb,
		input$file6use_cb,
		input$file7use_cb,
		input$file8use_cb);
	
	return(file_cb);
}

files_specified=function(input, output, session){

	file_path=c(
		ntb(input$file1_upld$datapath),
		ntb(input$file2_upld$datapath),
		ntb(input$file3_upld$datapath),
		ntb(input$file4_upld$datapath),
		ntb(input$file5_upld$datapath),
		ntb(input$file6_upld$datapath),
		ntb(input$file7_upld$datapath),
		ntb(input$file8_upld$datapath)
	);
	
	return(file_path);
}

function(input, output, session) {
	
	output$downloadNotReadyTxt=renderText({"Not Ready"});
	output$runNotReadyTxt=renderText({"Not Ready"});
	
	output$downloadReadyTxt=renderText({""});
	output$runReadyTxt=renderText({""});
	
	output$appstdoutReadyTxt=renderText({"(Nothing yet...)"});
	
	observeFileUploadsEvents(input, output, session);
	observeCheckboxEvents(input, output, session);

	observeEvent(input$runButton, {
		# This is where the uploaded file is on the server
		filelist=files_specified(input, output, session);
			
		# This is the name of the file the user uploaded (w/o the full path)
		filenames=c(
			ntb(input$file1_upld$name),
			ntb(input$file2_upld$name),
			ntb(input$file3_upld$name),
			ntb(input$file4_upld$name),
			ntb(input$file5_upld$name),
			ntb(input$file6_upld$name),
			ntb(input$file7_upld$name),
			ntb(input$file8_upld$name)
		);
		
		# This is where the uploaded file is on the server
		fileuse=files_selected(input, output, session);
		
		# Select out the files that can be used
		isNotEmpty=ifelse(filelist!="", T, F);
		filelist=filelist[fileuse&isNotEmpty];
		filenames=filenames[fileuse&isNotEmpty];
		
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
		output$appstdoutReadyTxt=renderText({paste(mydata, collapse="\n"); })
		
		output$downloadNotReadyTxt=renderText({""});
		output$downloadReadyTxt=renderText({"Ready"});
		
		session$userData[["output_files"]]=c(
			first=paste(output_root, ".accumulated.detail.tsv", sep=""),
			second=paste(output_root, ".accumulated.simplify.tsv", sep="")
		);
		
	});
	
	
	output$downloadButton=downloadHandler(
		filename=function(){
			paste("Results.zip");
		},
		content=function(filename){
			cat("Write to this file:\n");
			print(filename);
			
			# Zip file together and place results into the filename specified by shiny
			cat("Zipping these files:\n");
			print(session$userData[["output_files"]]);
			
			zip(zipfile=filename, files=session$userData[["output_files"]]);
			cat("Done zipping.\n");
		}
	);
	

}