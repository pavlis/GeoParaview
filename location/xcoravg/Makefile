BIN=xcoravg

#MAN1=		# man pages for commands
ldlibs=	$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)
OBJS=xcoravg.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
