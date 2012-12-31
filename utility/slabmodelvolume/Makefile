BIN=slabmodelvolume
#ldlibs=-lgeocoords -lcgeom -lgclgrid -lseispp $(DBLIBS) -lperf 
ldlibs=-lgeocoords -lcgeom -lgclgrid -lseispp $(TRLIBS) -lperf 

#ANTELOPEMAKELOCAL = $(ANTELOPE)/local/include/antelopemake.local
include $(ANTELOPEMAKE)  	
#include $(ANTELOPEMAKELOCAL)

OBJS=slabmodelvolume.o 

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
