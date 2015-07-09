BIN=pbpathtolatlon
ldlibs=-lgclgrid -lgeocoords -lseispp $(DBLIBS) -lperf 

include $(ANTELOPEMAKE)  	

OBJS=pbpathtolatlon.o  

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
