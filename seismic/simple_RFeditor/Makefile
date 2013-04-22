BIN=simple_RFeditor
cflags=-g
ldflags=-lseispp -lgclgrid -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE)
OBJS=	simple_RFeditor.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
