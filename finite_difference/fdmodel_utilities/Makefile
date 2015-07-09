BIN=sym_to_seispp
#ldlibs=-lseispp -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid
ldlibs=-lseispp -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lgclgrid -ltrvltm

#CXX=g++ -g
include $(ANTELOPEMAKE)  	
OBJS=sym_to_seispp.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
sym2gclgrid : sym2gclgrid.o
	$(CXX) $(CCFLAGS) -o $@ sym2gclgrid.o $(LDFLAGS) $(LDLIBS)
