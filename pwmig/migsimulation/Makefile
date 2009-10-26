
all clean Include install installMAN pf relink tags uninstall test :: FORCED
	@-if amakelocal XMOTIFLIB XMOTIFINCLUDE BOOSTLIB BOOSTINCLUDE; then \
	    $(MAKE) -f Makefile2 $@ ; \
	fi

FORCED:

