BIN=clone_cat
PF=clone_cat.pf

#MAN1=	clonse_cat.1
cflags=-g
ldflags= 	# flags to linker/loader
#  This list is longer than necessary because the get_required_pf function is loaded
# from libmultiwavelet.  I need a glp utility library for this type function
ldlibs=-lgenloc -lglputil -xlic_lib=sunperf,fui,fai,fsu $(DBLIBS) -ltrvltm -ldl -lm

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=clone_cat.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
