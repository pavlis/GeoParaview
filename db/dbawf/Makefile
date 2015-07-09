BIN=dbawf
PF=dbawf.pf
include $(ANTELOPEMAKE)
LDLIBS= -lscv2 -ldbl2 -lgrx $(DBLIBS) $(TRLIBS) -lX11  -lsunmath -lm
DIRS=
OBJS= dbawf.o read_wf.o trace_subs.o
CC=cc
$(BIN):$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) 
