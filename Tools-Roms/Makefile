include ../src/Makedefs.inc

ALL = partit ncjoin ncjoin_mpi join_child_bry nc3to4z zslice

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
	
partit: partit.o insert_node.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o lenstr.o $(LCDF)	

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

 ZSLICE = roms_read_write.o	zslice.o set_depth.o sigma_to_z_intr.o tiling.o\
          read_add_write.o read_string.o lenstr.o
          
zslice: $(ZSLICE)
	$(LDR) $(FFLAGS) -o zslice $(ZSLICE) $(LCDF)

clean:
	/bin/rm -f *.f *.o *.a *.s *.trace

clobber: clean
	/bin/rm -f $(ALL) *.bak
