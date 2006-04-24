BIN=uwitodb
PF=uwitodb.pf
MAN1=		# man pages for commands

cflags=-g
ldlibs=$(DBLIBS)

include $(ANTELOPEMAKE)
OBJS=	uwitodb.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
