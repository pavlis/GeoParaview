BIN=rpitomo2db
cflags=-g
ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) $(F77LIBS)

include $(ANTELOPEMAKE)
OBJS=rpitomo2db.o utmsubs.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
