LIB=libseispp.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
INCLUDE=metadata.h
MAN3=metadata.3
CC=CC -Kpic
AR=CC -xar

cflags=-g
ldlibs= -L/opt/SUNWspro/lib $(DBLIBS) -lperf

include $(ANTELOPEMAKE)  	
DLD=CC -G -Kpic


OBJS=metadata.o pfstream_interface.o
$(LIB) : $(OBJS)
	rm -f $@
#	$(AR) $(ARFLAGS) $@ $(LORDER) $(OBJS) $(TSORT)
	$(AR) -o $@ $(LORDER) $(OBJS) $(TSORT)
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(DLD) $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
	#CC -G -Kpic $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
