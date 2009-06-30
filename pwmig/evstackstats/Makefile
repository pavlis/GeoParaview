BIN=evstackstats

ldlibs=-lseispp $(DBLIBS)

include $(ANTELOPEMAKE)  	
OBJS=evstackstats.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
