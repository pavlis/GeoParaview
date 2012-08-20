BIN=db2rpitomo
#PF=rpitomo2db.pf
cflags=-g
ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) 

include $(ANTELOPEMAKE)
OBJS=db2rpitomo.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
