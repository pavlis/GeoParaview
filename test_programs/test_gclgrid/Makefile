BIN=testgcl
ldlibs=-lgclgrid -lseispp -lpfstream -lbrttutil -ltrvltm $(TRLIBS) $(DBLIBS) -lperf 

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 
OBJS=testgcl.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
