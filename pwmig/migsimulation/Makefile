BIN=migsimulation
PF=migsimulation.pf
ldlibs=-lseispp -lsynthetic -ltrvltm -lbrttutil $(TRLIBS) $(DBLIBS)  \
               -lperf -lgclgrid -lpwmig 

include $(ANTELOPEMAKE) 

OBJS=migsimulation.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
