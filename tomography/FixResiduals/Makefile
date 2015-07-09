BIN=fixresiduals
ldlibs=-lgclgrid -lseispp -lperf $(DBLIBS) 

include $(ANTELOPEMAKE)  	
OBJS=fixresiduals.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -g -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
