#!/usr/bin/Rscript


library(Bookmarker);
bm=bookmark_init("test.pdf");


bm=bookmark_insert(bm, "MAIN", "Bookmark Main", desc="This is an example of how to describe the main analysis document.");
bm=bookmark_insert(bm, "SECT", "First Section", pg_num=5, desc="This describes the first section");
bm=bookmark_insert(bm, "SECT", "Second Section", desc="This describes the second section");
bm=bookmark_insert(bm, "ITEM", "Item 1", pg_num=6, desc="Short description of item 1.");
bm=bookmark_insert(bm, "ITEM", "Item 2", pg_num=7, desc="Something quick about item 2.");
bm=bookmark_insert(bm, "SECT", "Third Section", desc="This describes the third section");

write_bookmarks_HTML(bm);
