BIN=sprt_band
PF=sprt_band.pf
CFLAGS = -g -I$(ANTELOPE)/include
F77LIBS=-L/opt/SUNWspro/lib -R/opt/SUNWspro/lib -lM77 -lF77 -lsunmath -lm
LDFLAGS = -L$(ANTELOPE)/lib 
ldlibs= -ldbl2 $(DBLIBS) $(TRLIBS) -lX11 
include $(ANTELOPEMAKE)
DBG=-g
OBJS= sprt_band.o load_station_table.o read_spectrum.o \
resample_spectrum.o
CC=cc
$(BIN):$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(F77LIBS)
