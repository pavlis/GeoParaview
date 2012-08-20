BIN=dbdefine_blasts

ldlibs=-lseispp -lgclgrid -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lm
include $(ANTELOPEMAKE)  	
OBJS=dbdefine_blasts.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
