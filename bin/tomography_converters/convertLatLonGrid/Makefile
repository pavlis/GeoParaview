BIN=convertLatLonGrid
cflags=-O3 
ldlibs=-lgclgrid -lseispp $(DBLIBS) -lperf

SUBDIR=/contrib
include $(ANTELOPEMAKE) 
include $(ANTELOPE)/contrib/include/antelopemake.local
CXXFLAGS += -I$(BOOSTINCLUDE)
LDFLAGS += -L$(BOOSTLIB)
OBJS=	convertLatLonGrid.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
