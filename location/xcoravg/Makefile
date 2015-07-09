BIN=xcoravg

#MAN1=		# man pages for commands
ldlibs=	-lseispp -ltrvltm -lgclgrid $(TRLIBS) $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE)
OBJS=xcoravg.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
