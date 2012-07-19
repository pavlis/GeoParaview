PF=gocad2vtk.pf
BIN=gocad2vtk
#ldlibs=-lgeocoords -lcgeom -lgclgrid -lseispp $(DBLIBS) -lperf 
#ldlibs=-lgeocoords -lcgeom -lgclgrid -lseispp $(TRLIBS) -lperf 
ldlibs=-lgeocoords -lgclgrid -lseispp $(TRLIBS) -lperf 

#ANTELOPEMAKELOCAL = $(ANTELOPE)/local/include/antelopemake.local
include $(ANTELOPEMAKE)  	
#include $(ANTELOPEMAKELOCAL)

OBJS=gocad2vtk.o 

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
