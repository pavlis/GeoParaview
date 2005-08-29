#makefile for MPI tutorial
#Bob Cruise - UITS-HPC - 08/23/01


##################################################

CC = mpicc
CFLAGS  = -g 
LDFLAGS = -lm

#CC= mpcc
#CC= mpcc_r

## We need to include VT.h for Vampirtrace profiling instrumentation. ##
## Compile with -g  and without optimization for TotalView debugger. ###
#CFLAGS = -g -I/libapps/Vampirtrace/include -qmaxmem=200000
#CFLAGS = -g #-bmaxdata:2000000000 
#CFLAGS = -g -b64 -q64


###### We need to link to the Vampirtrace library. #######
#LDFLAGS   = -lm -lld #-L/libapps/Vampirtrace/lib  -lVT
#LDFLAGS   = -lm -b64 -q64

PROG = elwave3d

#OBJS = basic_mpi.o getpar.o initialize.o model.o source.o abc.o surface.o communicate_mpi.o output.o update.o el3d_mpi.o 
OBJS = basic_mpi.o getpar.o read_setup.o initialize.o model.o source.o abc.o surface.o communicate_mpi.o output.o update.o el3d_mpi.o 
SRCS = $(OBJS:.o=.c)


all: $(PROG)

$(PROG):$(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f *.o
	rm -f core
