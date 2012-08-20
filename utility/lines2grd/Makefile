BIN=line2grd
#MAN1=line2grd.1
cflags=-g
ldlibs=$(DBLIBS) -lseispp -lgclgrid -lperf -lm -lseispp

include $(ANTELOPEMAKE)  	
OBJS=line2grd.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
