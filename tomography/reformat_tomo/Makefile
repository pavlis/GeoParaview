BIN=reformat_tomo
cflags=-g
ldlibs=-lsunmath -lm -lgclgrid $(F77LIBS)

include $(ANTELOPEMAKE)
OBJS=reformat_tomo.o utmsubs.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
