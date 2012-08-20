LIB=libhpcdb.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
BUNDLE=$(LIB:.a=.bundle) 
INCLUDE=ProcessingQueue.h DatascopeProcessingQueue.h HPCDBHandle.h
#MAN3=ProcessingQueue.3

include $(ANTELOPEMAKE)  
#
# The datascope stuff here should eventually end up in seispp
#
OBJS=DatascopeProcessingQueue.o HPCDBHandle.o 

$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  # the order "$(CXXARFLAGS) $@" is required
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(RM) $@
	$(CXXDLD) $(CXXDLDFLAGS) -o $(DLIB) $(OBJS) $(DLDLIBS) 

$(BUNDLE) : $(OBJS)
	$(RM) $@
	$(DLD)  -flat_namespace -bundle -undefined suppress -o $@   $(OBJS) $(LDFLAGS) $(DLDLIBS)
