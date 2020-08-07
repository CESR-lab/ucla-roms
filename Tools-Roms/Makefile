#include Makedefs.GCC
include Makedefs.IntelEM64T

ALL = cross_matrix chkindent mreplace partit ncjoin

all: $(ALL)

install: $(ALL)
	mv $(ALL) ~/bin

# Note: the following .f: dependencies are needed to make everything
# go through CPP and take custom defined CPPFLAGS rather than compiling
# executables from .F directly.

cross_matrix: cross_matrix.f
chkindent: chkindent.f

# The following two are built from the same source
# and differ only CPP macro setting:

mreplace.f: mreplace.F
	$(CPP) -P $(CPPFLAGS) -DREPLACE mreplace.F > mreplace.f
	

# Need to set up partit properly - previously in main roms makefile	
#partit: partit.o insert_node.o lenstr.o
#	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o lenstr.o $(LCDF)	

# Need to set up ncjoin properly - previously in main roms makefile
#ncjoin: ncjoin.o lenstr.o
#	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin ncjoin.o lenstr.o $(LCDF)	

depend: ./cross_matrix
	./cross_matrix *.F


clean:
	/bin/rm -f *.f *.o *.a *.s *.trace *~ $(COMP_FILES)

clobber: clean
	/bin/rm -f $(ALL) *.bak

# Need to set up Make.depend properly. cross_matrix didn't work for some reason.
# Make.depend was created but contained no dependencies.
include Make.depend
