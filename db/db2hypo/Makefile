BIN=db2artm
#MAN1=
ldlibs=-lseispp -ltrvltm $(DBLIBS) -lperf -lm
include $(ANTELOPEMAKE)
OBJS=db2artm.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
