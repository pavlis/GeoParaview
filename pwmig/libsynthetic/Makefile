LIB=libsynthetic.a
DLIB=$(LIB:.a=$(DSUFFIX))
BUNDLE=$(LIB:.a=.bundle)
INCLUDE=SyntheticSeismogram.h PointSourcePSSynthetic.h \
        vectorcls.h SimplePSPrimarySynthetic.h StaVariableLayeredSynthetic.h \
        SphericalRayPathArray.h   
cflags=-g
ldlibs=$(F77LIBS)

include $(ANTELOPEMAKE)
OBJS=	PointSourcePSSynthetic.o SimplePSPrimarySynthetic.o \
        SphericalRayPathArray.o StaVariableLayeredSynthetic.o \
        vectorcls.o


$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(RM) $@
	$(CXXDLD) $(CXXDLDFLAGS) -o $(DLIB) $(OBJS) $(DLDLIBS) 

$(BUNDLE) : $(OBJS)
	$(RM) $@
	$(DLD)  -flat_namespace -bundle -undefined suppress -o $@   $(OBJS) $(LDFLAGS) $(DLDLIBS)
