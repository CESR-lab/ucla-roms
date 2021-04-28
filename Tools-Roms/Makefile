include ../src/Makedefs.inc

ALL = partit ncjoin ncjoin_mpi join_child_bry

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
	
ncjoin_mpi: ncjoin_mod.o ncjoin_mpi.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin_mpi ncjoin_mod.o ncjoin_mpi.o lenstr.o $(LCDF)
	
nc3to4z: nc3to4z.o lenstr.o read_string.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o nc3to4z nc3to4z.o lenstr.o read_string.o $(LCDF)
	
mpc: mpc.F Makefile
	$(CPP) $(CPPFLAGS) -P mpc.F > ./mpc.f
	$(CFT) $(FFLAGS) $(LDFLAGS) -o mpc ./mpc.f
	
join_child_bry: ncjoin_mod.o join_child_bry.o  lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o join_child_bry ncjoin_mod.o join_child_bry.o lenstr.o  $(LCDF)

clean:
	/bin/rm -f *.f *.o *.a *.s *.trace

clobber: clean
	/bin/rm -f $(ALL) *.bak
