BIN=pwstack
PF=pwstack.pf
DATADIR=schemas
DATA=pwmig1.1

#MAN1=pwstack.1
cflags=-g
ldflags=-lseispp -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid

include $(ANTELOPEMAKE) 
OBJS=pwstack.o \
	apsubs.o \
	aperture.o \
	out.o \
	pwstack_process.o 
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
