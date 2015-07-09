BIN=slabproject
PF=slabproject.pf
#MAN1=pathtosurface.1
cflags=-g
ldlibs=$(DBLIBS) -lseispp -lgclgrid -lperf -lm -lseispp

include $(ANTELOPEMAKE)  	
OBJS=slabproject.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
