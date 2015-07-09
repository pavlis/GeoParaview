PF=cpsm.pf
BIN=cpsm
ldlibs=-lgeocoords -lcgeom -lgclgrid -lseispp $(DBLIBS) -lperf 

#ANTELOPEMAKELOCAL = $(ANTELOPE)/local/include/antelopemake.local
include $(ANTELOPEMAKE)  	
#include $(ANTELOPEMAKELOCAL)

OBJS=cpsm.o 

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
