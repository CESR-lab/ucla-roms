include ../src/Makedefs.inc

ALL = partit ncjoin

all:	$(info MAKING ROMS TOOLS:) mpc $(ALL)

#install: $(ALL)
#	mv $(ALL) ~/bin

# Note: the following .f: dependencies are needed to make everything
# go through CPP and take custom defined CPPFLAGS rather than compiling
# executables from .F directly.

chkindent: chkindent.f

# The following two are built from the same source
# and differ only CPP macro setting:

mreplace.f: mreplace.F
	$(CPP) -P $(CPPFLAGS) -DREPLACE mreplace.F > mreplace.f
	

# Need to set up partit properly - previously in main roms makefile	
partit: partit.o insert_node.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o lenstr.o $(LCDF)	

# Need to set up ncjoin properly - previously in main roms makefile
ncjoin: ncjoin.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin ncjoin.o lenstr.o $(LCDF)	
	
mpc: mpc.F Makefile
	$(CPP) $(CPPFLAGS) -P mpc.F > ./mpc.f
	$(CFT) $(FFLAGS) $(LDFLAGS) -o mpc ./mpc.f

clean:
	/bin/rm -f *.f *.o *.a *.s *.trace *~ $(COMP_FILES)

clobber: clean
	/bin/rm -f $(ALL) *.bak
