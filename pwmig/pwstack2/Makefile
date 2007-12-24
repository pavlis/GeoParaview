BIN=pwstack
PF=pwstack.pf
#DATADIR=schemas
#DATA=pwmig1.2

#MAN1=pwstack.1
cflags=-g
# Use the following for Solaris
ldflags=-lseispp -lgclgrid -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid -lseispp
# Following is required to work with g++ 3.2.3
#ldflags=-shared -lseispp -lgclgrid -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid


include $(ANTELOPEMAKE) 
OBJS=pwstack.o \
	PwmigFileHandle.o \
	apsubs.o \
	aperture.o \
	pwstack_process.o \
	compute_stack_coherence.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(DBX) $(LDFLAGS) $(LDLIBS)
