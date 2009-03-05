BIN=telecluster
PF=telecluster.pf

MAN1=telecluster.1
ldlibs= -lseispp -ltrvltm -lgclgrid $(TRLIBS) $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE) 
OBJS= telecluster.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
