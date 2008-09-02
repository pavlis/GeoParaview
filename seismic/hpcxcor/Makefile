BIN=hpcxcor
#PF=hpcxcor.pf
#MAN1=hpcxcor.1

ldlibs=-lseispp -ltrvltm -lgclgrid $(TRLIBS)  $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE) 

# all but hpcxcor need to eventually be removed. AnalysisSetting and XcorProcessingEngine
# need to be moved to seispp
OBJS=hpcxcor.o 

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(cxxflags) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
