BIN=discard_arrivals
#PF=discard_arrivals.pf

ldflags= 	# flags to linker/loader
#  This list is longer than necessary because the get_required_pf function is loaded
# from libmultiwavelet.  I need a glp utility library for this type function
ldlibs= $(DBLIBS) 

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=discard_arrivals.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
