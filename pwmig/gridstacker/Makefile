BIN=gridstacker
DIRS=
PF=gridstacker.pf

#MAN1=gridstacker.1
ldflags=-lseispp -lgclgrid $(TRLIBS) $(DBLIBS) -lperf -lseispp -lm

include $(ANTELOPEMAKE) 
OBJS= gridstacker.o \
  GridScratchFileHandle.o \
  GridStackPenaltyFunction.o
mdutils.cc : pwmig.h
pwmig.cc : pwmig.h
ray_transformation_operator.cc : pwmig.h
raygrid.cc : pwmig.h
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
