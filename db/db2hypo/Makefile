BIN=db2artm
#MAN1=
fflags=$(DBLIBS)
include $(ANTELOPEMAKE)
OBJS=db2artm.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
