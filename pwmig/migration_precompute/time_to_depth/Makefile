BIN=time_to_depth
CSSFLAGS= -I/opt/antelope/4.5/include

cflags=-g 
ldlibs= -L/opt/SUNWspro/lib $(DBLIBS) 

include $(ANTELOPEMAKE)
CC=CC
CXX=CC
LD=CC
OBJS=time_to_depth.o vmod1d.o raytrace.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
