BIN=vtk_gcl_converter

ldlibs=-lm -lgclgrid -lseispp -lperf $(DBLIBS) 

include $(ANTELOPEMAKE)
OBJS=vtk_gcl_converter.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
