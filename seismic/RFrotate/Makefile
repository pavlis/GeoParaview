BIN=RFrotate
MAN1=RFrotate.1

ldlibs=-lseispp -lgclgrid -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm
include $(ANTELOPEMAKE)
OBJS=RFrotate.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
