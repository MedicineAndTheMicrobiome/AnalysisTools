
bookmark_init=function(assoc_content_name){
	bm=list()
	
	bm[["content_name"]]=assoc_content_name;
	bm[["last_page"]]=0;
	bm[["bookmarks"]]=list();
	bm[["bookmark_idx"]]=0;
	
	return(bm);
}

bookmark_insert=function(bm, level="ITEM", name, pg_num=-1, desc=""){
	# Used to describe and/or introduce the entire
	# analysis.  There should only be one, but you
	# could put in multiple, if desired.

	if(!any(level==c("MAIN", "SECT", "ITEM"))){
		cat("Error:  Bookmark level type must be one of MAIN, SECT, or ITEM\n");
		quit(-1);
	}

	if(pg_num==""){
		pg_num=bm[["_lastpage"]]+1;
	}


	bmix=bm[["bookmark_idx"]]+1;
	bm[["bookmarks"]][[bmix]]=c(name, desc, pg_num, level);
	bm[["bookmark_idx"]]=bmix;

	bm[["_lastpage"]]=pg_num;

	return(bm);
}

get_last_page=function(bm){
	return(bm[["_lastpage"]]);
}

write_bookmarks_HTML=function(bm){
	
	bm_fname=paste(bm[["content_name"]], ".html", sep="");
	fh=file(bm_fname, "w");

	cat(file=fh, "<html>\n");
	cat(file=fh, "<head></head>\n");
	cat(file=fh, "<body>\n");

	link=function(pgnum, link_text){
		if(pgnum==-1){
			# If there's no page number, then it's unassociated text
			return(link_text);
		}
		return(
			paste("<a href=\"", bm[["content_name"]], "#page=", pgnum, "\" target=\"frame2\">",
				link_text, "</a>", sep="")

		);
	}

	last_bm_ix=bm[["bookmark_idx"]];
	if(last_bm_ix==0){
		cat("No bookmarks made.\n");
		cat(file=fh, "No bookmarks...\n");
	}else{
		for(i in 1:last_bm_ix){
			bm_entry=bm[["bookmarks"]][[i]];

			cat("Working on:\n");
			print(bm_entry);

			if(bm_entry[4]=="MAIN"){
				str=paste("<h1>", link(bm_entry[3], bm_entry[1]), "</h1>\n<p>", bm_entry[2], "</p>", sep="");

			}else if(bm_entry[4]=="SECT"){
				str=paste("<h2>", link(bm_entry[3], bm_entry[1]), "</h2>\n<p>", bm_entry[2], "</p>", sep="");

			}else if(bm_entry[4]=="ITEM"){
				str=paste("<li>", link(bm_entry[3], bm_entry[1]), ": ", bm_entry[2], "</li>");
			}

			cat(file=fh, str, "\n", sep="");
		}
	}


	cat(file=fh, "</body>\n");
	cat(file=fh, "</html>\n");

	cat("HTML written.\n");
}





