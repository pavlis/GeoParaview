BIN=rpitomo2dx
cflags=-g
ldlibs=-lsunmath -lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) $(F77LIBS)

include $(ANTELOPEMAKE)
OBJS=rpitomo2dx.o utmsubs.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
