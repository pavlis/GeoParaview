BIN=flowtotimelines
#MAN1=flowtotimelines.1
cflags=-g
ldlibs=-lgeocoords -lgclgrid -lperf -lm -lseispp

include $(ANTELOPEMAKE)  	
OBJS=flowtotimelines.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
