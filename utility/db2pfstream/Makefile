BIN=db2pfstream
LIB=
DLIB=$(LIB:.a=$(DSUFFIX))
DATADIR=
DATA=
INCLUDE=
PF=

MAN1=db2pfstream.1
MAN3=
MAN5=

CLEAN=

cflags=
fflags=
ldflags=
ldlibs=$(DBLIBS)

include $(ANTELOPEMAKE)

DIRS=



OBJS=db2pfstream.o pfget_Attribute_map.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
