BIN=plot_spectra
PF=plot_spectra.pf
ldlibs = -lmap -lgrx $(DBLIBS) -lX11 
include $(ANTELOPEMAKE)
#F77LIBS=-L/opt/SUNWspro/lib -R/opt/SUNWspro/lib -lM77 -lF77 
F77LIBS= -lM77 -lF77 -L/opt/SUNWspro/lib -lm
DBG=-g
OBJS= plot_spectra.o read_spectrum.o init_plot.o llaxis.o
CC=cc
$(BIN):$(OBJS)
	f77 $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(F77LIBS)
