BIN=surf_vtk_converter
ldlibs=-lseispp -ltrvltm -lgclgrid -lseispp $(TRLIBS) $(DBLIBS) -lperf -lwffil -lbrttutil -lm
SUBDIR=/contrib

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)         
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I/opt/boost/include/
OBJS= surf_vtk_converter.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
