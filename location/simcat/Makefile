BIN=simcat
PF=simcat.pf

#MAN1=	simcat.1
ldflags= 	# flags to linker/loader
#  This list is longer than necessary because the get_required_pf function is loaded
# from libmultiwavelet.  I need a glp utility library for this type function
ldlibs=-lgenloc -lgclgrid -lglputil -xlic_lib=sunperf -lglputil $(DBLIBS) -ltt -ldl -lm

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=simcat.o load_stations.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
