BIN=syndecon
CC=CC
CXX=CC
LD=CC
CSSFLAGS= -I/opt/antelope/4.4/include

cflags=-g 
ldlibs= -L/opt/SUNWspro/lib -lseispp -lpfstream -lbrttutil $(DBLIBS) -lperf 
include $(ANTELOPEMAKE)
DIRS=
OBJS=syndecon.o sac.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
