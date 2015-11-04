BIN=geo2cart

SUBDIR=/contrib
ldlibs=-lm -lseispp -lgclgrid -lperf $(DBLIBS) 

include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
CXXFLAGS += -g
CXXFLAGS += -I/opt/boost/include/
OBJS=geo2cart.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
