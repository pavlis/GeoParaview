LIB=libseispp.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
INCLUDE=metadata.h seispp.h dmatrix.h
MAN3=metadata.3 dmatrix.3

cflags=-g
ldlibs= -L/opt/SUNWspro/lib $(DBLIBS) -lperf

include $(ANTELOPEMAKE)  	
DLD=CC -G -Kpic
CC=CC -Kpic
AR=CC -xar


OBJS=metadata.o create_destroy.o dmatrix.o operators.o pfstream_interface.o \
rotate_procedures.o rotation.o
$(LIB) : $(OBJS)
	rm -f $@
#	$(AR) $(ARFLAGS) $@ $(LORDER) $(OBJS) $(TSORT)
	$(AR) -o $@ $(LORDER) $(OBJS) $(TSORT)
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(DLD) $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
	#CC -G -Kpic $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
