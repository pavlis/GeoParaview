BIN=pwstack
PF=pwstack.pf
CC=CC
LD=CC

#MAN1=pwstack.1
cflags=-g
ldflags=$(DBLIBS)

include $(ANTELOPEMAKE) 
OBJS=pwstack.o \
	apsubs.o \
	aperture.o \
	create_destroy.o \
	pwstack_process.o 
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
