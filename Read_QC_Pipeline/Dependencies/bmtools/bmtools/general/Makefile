CC = gcc
CXX = g++
DEBUG = 
FLAGS = $(DEBUG) -O6
CXXFLAGS = $(FLAGS) -DSTANDALONE
CFLAGS = $(FLAGS) -DSTANDALONE
LDFLAGS = $(FLAGS)
LIBS =

ALL_SRC = $(libof_general_SRC)
ALL_OBJ = $(ALL_SRC:%=%.o)
ALL_DEP = $(ALL_SRC:%=%.d)

libof_general_SRC = \
	cseqcoding \
	cpermutator8b \
	cprogressindicator \
	cshortreader \
	csamrecords \
	capp

libof_general_OBJ = $(libof_general_SRC:%=%.o)

%.o: %.cpp %.d
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.d: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ -MM $<


all: libof-general.a

test:

clean:
	-rm $(ALL_OBJ) $(ALL_DEP) libof-general.a

depend: $(ALL_DEP)

libof-general.a: $(libof_general_OBJ)
	$(AR) r $@ $(libof_general_OBJ)

-include $(ALL_DEP)
