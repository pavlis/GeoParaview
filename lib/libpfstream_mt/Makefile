LIB=libpfstream.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
DATADIR=
DATA=
INCLUDE=pfstream.h
MAN3=
MAN5=		# man pages for file formats

cflags=-g
ldlibs= -lbrttutil $(THREAD) $(DBLIBS) 

include $(ANTELOPEMAKE)  	


OBJS=pfstream.o pfdup.o pfget_Attribute_map.o pthread_routines.o
$(LIB) : $(OBJS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(LORDER) $(OBJS) $(TSORT)
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(DLD) $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
