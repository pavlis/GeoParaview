BIN=pmelavg
PF=pmelavg.pf

#MAN1=	pmelavg.1
ldflags=-L/opt/SUNWspro/lib
ldlibs= -lgenloc -lglputil $(DBLIBS) $(F77LIBS) \
-ltrvltm -ldl -lm -xlic_lib=sunperf,fui,fai,fsu -dalign
include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=pmelavg.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
