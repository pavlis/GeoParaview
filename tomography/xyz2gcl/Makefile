BIN=xyz2gcl
PF=xyz2gcl.pf
cflags=-g
ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) 

include $(ANTELOPEMAKE)
OBJS=xyz2gcl.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
