BIN=db_ratio
cflags = -g
ldlibs = -lmap -lgrx $(DBLIBS) -lX11 -lm
F77LIBS=-L/opt/SUNWspro/lib -R/opt/SUNWspro/lib -lM77 -lF77 -lsunmath -lm
include $(ANTELOPEMAKE)
DBG=-g
OBJS= db_ratio.o data_gather.o load_station_table.o median_float.o read_spectrum.o save_spectrum.o init_plot.o llaxis.o
CC=cc
$(BIN):$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(F77LIBS)
