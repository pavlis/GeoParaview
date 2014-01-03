BIN=volume_counter
PF=volume_counter.pf
MAN1=volume_counter.1
cflags=-g
ldlibs=$(DBLIBS) -lseispp -lgclgrid -lperf -lm -lseispp

include $(ANTELOPEMAKE)  	
OBJS=volume_counter.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
