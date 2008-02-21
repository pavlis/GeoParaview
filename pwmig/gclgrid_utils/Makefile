BIN=gridsum
#MAN1=pwmig.1
cflags=-g
ldlibs=-lseispp -ltrvltm -lgclgrid $(TRLIBS) $(DBLIBS) -lperf -lm -lseispp

include $(ANTELOPEMAKE) 
gridsum : gridsum.o
	$(CXX) $(CXXFLAGS) -o gridsum gridsum.o $(LDFLAGS) $(LDLIBS)
dumpgridcoords : dumpgridcoords.o
	$(CXX) $(CXXFLAGS) -o dumpgridcoords  dumpgridcoords.o  $(LDFLAGS) $(LDLIBS)
