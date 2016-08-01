CC = gcc
CXX = g++
DEBUG = 
FLAGS = $(DEBUG) -O6
CXXFLAGS = $(FLAGS) -DSTANDALONE -I ../general
CFLAGS = $(FLAGS) -DSTANDALONE -I ../general
LDFLAGS = $(FLAGS) -L ../general
LIBS = $(general_LIBS:%=-l%) $(ADD_LIBS:%=-l%)

ADD_LIBS=

ifdef SRA_VERSION
VDB_INCLUDES ?= ERROR
FLAGS += -DWITH_SRA -I$(VDB_INCLUDES) 
LDFLAGS += -L$(VDB_LIBRARIES)
ADD_LIBS += \
    srareader \
	sradb \
	srapath \
	vdb \
	kdb \
	vfs \
	krypto \
	kfg \
    kfs \
    klib \
    kproc \
    kq \
    ksproc \
    ksrch \
    kurl \
    kxml \
    m z bz2 pthread
endif

general_LIBS = of-general

ALL_SRC = $(bmtool_SRC) $(bmfilter_SRC) $(bmdiff_SRC) $(bmdump_SRC)
ALL_OBJ = $(ALL_SRC:%=%.o)
ALL_DEP = $(ALL_SRC:%=%.d)

bmtool_SRC = bmtool-main
	
bmtool_OBJ = $(bmtool_SRC:%=%.o)

bmdiff_SRC = bmdiff-main
	
bmdiff_OBJ = $(bmdiff_SRC:%=%.o)

bmdump_SRC = bmdump-main
	
bmdump_OBJ = $(bmdump_SRC:%=%.o)

bmfilter_SRC = \
	creadtagger \
	bmfilter-main
	
bmfilter_OBJ = $(bmfilter_SRC:%=%.o)

%.o: %.cpp %.d
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.d: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ -MM $<


all: bmtool bmfilter bmdiff bmdump extract_fullseq

test: bmtool bmfilter bmdiff bmdump
	./bmtool -hV
	./bmfilter -hV
	./bmdiff -hV
	./bmdump -hV

clean:
	-rm $(ALL_OBJ) $(ALL_DEP) bmtool bmfilter bmdiff bmdump extract_fullseq

depend: $(ALL_DEP)

bmtool: $(bmtool_OBJ) Makefile $(general_LIBS:%=../general/lib%.a)
	$(CXX) -o $@ $(bmtool_OBJ) $(LDFLAGS) $(LIBS)

bmdiff: $(bmdiff_OBJ) Makefile $(general_LIBS:%=../general/lib%.a)
	$(CXX) -o $@ $(bmdiff_OBJ) $(LDFLAGS) $(LIBS)

bmdump: $(bmdump_OBJ) Makefile $(general_LIBS:%=../general/lib%.a)
	$(CXX) -o $@ $(bmdump_OBJ) $(LDFLAGS) $(LIBS)

bmfilter: $(bmfilter_OBJ) Makefile $(general_LIBS:%=../general/lib%.a)
	$(CXX) -o $@ $(bmfilter_OBJ) $(LDFLAGS) $(LIBS)

extract_fullseq: extract_fullseq.c Makefile
	$(CC) -o $@ extract_fullseq.c

-include $(ALL_DEP)
