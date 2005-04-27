BIN=gclfield2vtk

vtkroot =/space/pavlis/VTK
vtkcore =-L$(vtkroot)/lib/vtk -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon
vtkext =  -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtkftgl -lvtkfreetype
vtkgl =  -lGL
vtklibs = $(vtkcore) $(vtkext) $(vtkgl)
vtkincs=-I$(vtkroot)/include/vtk

cflags=-g $(vtkincs)
ldlibs=-lsunmath -lm -lgclgrid -lseispp -lperf $(DBLIBS) $(TRLIBS) $(F77LIBS) $(vtklibs)

include $(ANTELOPEMAKE)
OBJS=gclfield2vtk.o vtk_output.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
