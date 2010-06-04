BIN=pbpathtolatlon
ldlibs=-lgclgrid -lseispp $(DBLIBS) -lperf 

include $(ANTELOPEMAKE)  	

OBJS=pbpathtolatlon.o   PlateBoundaryPath.o 

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
