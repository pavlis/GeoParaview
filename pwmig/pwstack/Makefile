BIN=pwstack
PF=pwstack.pf

#MAN1=pwstack.1
cflags=-g
ldflags=-lglputil -lpfstream -lseispp $(DBLIBS) -lgclgrid -lperf -lsunmath -lm

include $(ANTELOPEMAKE) 
OBJS=pwstack.o \
	apsubs.o \
	aperture.o \
	pwstack_process.o 
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
