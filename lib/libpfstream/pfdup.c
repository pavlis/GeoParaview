
#include <stdio.h>
#include <stdarg.h>
#include "stock.h"


Pf             *
pfdup (Pf * old)
{
    Pf             *new,
                   *value;
    char           *key;
    Tbl            *tbl;
    int             i,
                    n;

    new = pfnew (old->type);
    switch (old->type) {
      case PFSTRING:
	new->value.s = strdup ((char *) old->value.s);
	break;

      case PFFCT:
	new->value.pffunct->s = strdup ((char *) old->value.s);
	break;

      case PFARR:
      case PFFILE:
	tbl = keysarr (old->value.arr);
	n = maxtbl (tbl);
	for (i = 0; i < n; i++) {
	    key = (char *) gettbl (tbl, i);
	    value = (Pf *) getarr (old->value.arr, key);
	    setarr (new->value.arr, key, pfdup (value));
	}
	freetbl (tbl, 0);
	break;

      case PFTBL:
	tbl = old->value.tbl;
	n = maxtbl (tbl);
	for (i = 0; i < n; i++) {
	    value = (Pf *) gettbl (tbl, i);
	    pushtbl (new->value.tbl, pfdup (value));
	}
	break;

      default:
	die (0, "Bad type in pfdup : %d\n", old->type);
	break;
    }
    return new;
}
