LIB=libpfstream.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
DATADIR=
DATA=
INCLUDE=pfstream.h
MAN3=
MAN5=		# man pages for file formats

cflags=-g
ldlibs= $(DBLIBS) 

include $(ANTELOPEMAKE)  	


OBJS=pfstream.o pfdup.o pfget_Attribute_map.o
$(LIB) : $(OBJS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(LORDER) $(OBJS) $(TSORT)
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(DLD) $(DLDFLAGS) $(LDLIBS) -o $(DLIB) $(LORDER) $(OBJS) $(TSORT)
