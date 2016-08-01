#SUBDIRS = general oalign tests oligofar bmtagger tex man
SUBDIRS = general bmtagger



all depend test clean:
	set -e ; for file in $(SUBDIRS) ; do make -C $$file $(MAKEFLAGS:%=-%) $@ ; done


