BIN=computeresiduals
#MAN1=computeresiduals.o
ldlibs=-lgclgrid -lseispp -lperf $(DBLIBS) 

include $(ANTELOPEMAKE)  	
OBJS=computeresiduals.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -g -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
