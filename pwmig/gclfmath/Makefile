BIN=gclfmath

ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) 

include $(ANTELOPEMAKE)
OBJS=gclfmath.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
