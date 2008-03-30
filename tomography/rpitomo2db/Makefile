BIN=rpitomo2db
PF=rpitomo2db.pf
cflags=-g
ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) 

include $(ANTELOPEMAKE)
OBJS=rpitomo2db.o utmsubs.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
