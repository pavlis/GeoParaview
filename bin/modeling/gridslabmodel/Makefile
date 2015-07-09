BIN=gridslabmodel


ldlibs= $(DBLIBS) -lgclgrid -lseispp -ltrvltm $(TRLIBS) -lperf

include $(ANTELOPEMAKE)

OBJS=gridslabmodel.o
$(BIN) : $(OBJS)
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
