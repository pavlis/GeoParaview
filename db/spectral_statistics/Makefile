BIN=spectral_statistics
PF=spectral_statistics.pf
CFLAGS = -g -I$(ANTELOPE)/include
F77LIBS=-L/opt/SUNWspro/lib -R/opt/SUNWspro/lib -lM77 -lF77 -lsunmath -lm
LDFLAGS = -L$(ANTELOPE)/lib 
ldlibs= -ldbl2 $(DBLIBS) $(TRLIBS) -lX11 
include $(ANTELOPEMAKE)
DBG=-g
OBJS= spectral_statistics.o calc_statistics.o read_spectrum.o \
resample_spectrum.o
CC=cc
$(BIN):$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(F77LIBS)
