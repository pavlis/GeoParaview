BIN=dbmwspec
DATADIR=schemas2
DATA=spectrum1.0
PF=dbmwspec.pf
MAN1=dbmwspec.1
cflags=	-g
fflags= 
ldflags= 
ldlibs= -lsunmath -lmw -lscv2 -ldbl2 -lgrx $(DBLIBS) $(F77LIBS) $(TRLIBS) -lX11 -lm
include $(ANTELOPEMAKE)  

OBJS=dbmwspec.o \
	correct_for_response.o \
	powspc.o \
	save_spectrum.o \
	trace_subs.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
