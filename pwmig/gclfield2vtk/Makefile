BIN=gclfield2vtk
PF=gclfield2vtk.pf

vtkroot =/space/pavlis/VTK
vtkcore =-L$(vtkroot)/lib/vtk -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon
vtkext =  -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtkftgl -lvtkfreetype
vtkgl =  -lGL
vtklibs = $(vtkcore) $(vtkext) $(vtkgl)
vtkincs=-I$(vtkroot)/include/vtk

cflags=-g $(vtkincs)
ldlibs=-lsunmath -lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) $(F77LIBS) $(vtklibs)
#ldlibs=-lsunmath -lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) $(vtklibs)

#include antelopemake.erde
include $(ANTELOPEMAKE)
OBJS=gclfield2vtk.o vtk_output.o vtk_output_GCLgrid.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
