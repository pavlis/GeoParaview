BIN=convertDNA13
cflags=-g 
ldlibs=-lgclgrid -lseispp $(DBLIBS) -lperf

include $(ANTELOPEMAKE) 
OBJS=	convertDNA13.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
