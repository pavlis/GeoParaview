BIN=extract_events
PF=extract_events.pf
DBG=-g
ldlibs= -lgclgrid -lseispp $(TRLIBS) -lperf

include $(ANTELOPEMAKE)
OBJS=	extract_events.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
