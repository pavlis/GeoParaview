BIN=convertchen
cflags=-g 
ldlibs=-lgclgrid -lseispp $(DBLIBS) -lperf

include $(ANTELOPEMAKE) 
OBJS=	convertchen.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
