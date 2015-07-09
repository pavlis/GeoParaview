DATADIR=schemas
DATA=xcor1.0
BIN=xcor_proto
ldlibs=$(X11LIBS) -lseispp -lseisppplot $(TRLIBS) -lXm -lXt -lX11 -lm

include $(ANTELOPEMAKE)  	
OBJS=xcor_proto.o AnalysisSetting.o XcorProcessingEngine.o
XCOR=AnalysisSetting.o XcorProcessingEngine.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
xcor_analyzer : xcor_analyzer.o $(XCOR)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ xcor_analyzer.o $(XCOR) $(LDFLAGS) $(LDLIBS)
xcor_preprocessor : xcor_preprocessor.o $(XCOR)
	$(RM) $@
#	purify  $(CXX) $(CCFLAGS) -o $@ xcor_preprocessor.o $(XCOR) $(LDFLAGS) $(LDLIBS)
	$(CXX) $(CCFLAGS) -o $@ xcor_preprocessor.o $(XCOR) $(LDFLAGS) $(LDLIBS)
allxcor : $(BIN) xcor_analyzer xcor_preprocessor
