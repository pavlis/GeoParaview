BIN=migsimulation
PF=migsimulation.pf
#MAN1=migsimulation.1
ldlibs=-lseispp -ltrvltm -lbrttutil $(TRLIBS) $(DBLIBS)  \
               -lperf -lgclgrid -lpwmig 

include $(ANTELOPEMAKE) 

OBJS=migsimulation.o SimplePSPrimarySynthetic.o PointSourcePSSynthetic.o SphericalRayPathArray.o vectorcls.o
#OBJS=migsimulation.o SimplePSPrimarySynthetic.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

#clean:
#	-rm ${BIN}
#	-rm ${OBJS}
