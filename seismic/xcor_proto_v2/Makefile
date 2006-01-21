DATADIR=schemas
DATA=xcor1.0
BIN=xcor_proto
ldlibs=$(X11LIBS) -lseispp -lseisppplot $(TRLIBS) -lXm -lXt -lX11 -lm

include $(ANTELOPEMAKE)  	
OBJS=xcor_proto.o AnalysisSetting.o XcorProcessingEngine.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
