BIN=pathtosurface
PF=pathtosurface.pf
#MAN1=pathtosurface.1
cflags=-g
ldlibs=$(DBLIBS) -lseispp -lgclgrid -lperf -lm -lseispp

include $(ANTELOPEMAKE)  	
OBJS=pathtosurface.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
