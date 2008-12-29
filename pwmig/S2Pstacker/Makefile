BIN=S2Pstacker

DATADIR=schemas/css3.0.ext
DATA=pwstack

cflags=-g
# Use the following for Solaris
ldlibs= -lseispp -lgclgrid -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid -lseispp
# Following is required to work with g++ 3.2.3
#ldflags=-shared -lseispp -lgclgrid -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid


include $(ANTELOPEMAKE) 
OBJS=S2Pstacker.o \
	apsubs.o \
	aperture.o \
	pwstack_process.o \
	compute_stack_coherence.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
