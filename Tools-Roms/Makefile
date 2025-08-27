# Makefile for UCLA ROMS tools
#========== ======= =========== ======== === ==== =====
# This makefile requires that you have the following 
# environment variables defined: 
#  ROMS_ROOT, NETCDF_ROOT, NETCDFC_ROOT, MPIHOME, FCOMP
#
#  Example entries in a .bashrc file:
#
# export FCOMP=mpiifx
# export MPIHOME=/opt/intel/oneapi/mpi/latest
# export NETCDF_ROOT=/opt/netcdf-fortran-4.6.1
# export NETCDFC_ROOT=/opt/netcdf-c-4.9.2
#
# export LD_LIBRARY_PATH=${MPIHOME}/lib/:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=${NETCDFC_ROOT}/lib/:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=${NETCDF_ROOT}/lib/:$LD_LIBRARY_PATH

# export ROMS_ROOT=/home/username/ucla-roms
# export PATH=${ROMS_ROOT}/Tools-Roms:$PATH
# export PATH=${NETCDFC_ROOT}/bin:$PATH
# export PATH=${NETCDF_ROOT}/bin:$PATH

#export MARBL_ROOT=/home/username/MARBL
#export USEMARBL=TRUE (if the MARBL library is available)

#
.SUFFIXES: .o .f .F .f90


.F.o:
	@$(CPP) -P $*.F | mpc > $*.f
	$(FCOMP) -c $(FFLAGS) $(ICDF) -o $*.o $*.f
	@/bin/rm -f $*.f
.F.f:
	$(CPP) -P  $*.F | mpc > $*.f
.f.o:
	$(FCOMP) -c $(FFLAGS) $(ICDF) -o $*.o $*.f
.f90.o:
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f90

     CPP = cpp -traditional

# Path names of NetCDF and MPI libraries:

	NETCDF_INC = $(NETCDF_ROOT)/include
	NETCDF_LIB = -L$(NETCDF_ROOT)/lib -lnetcdff

ifeq ($(FCOMP),mpiifx)
	MPI_INC = ${MPIHOME}/include/mpi
else
	MPI_INC = ${MPIHOME}/include
endif


ifeq ($(FCOMP),mpiifx)
	CPATH = ${MPI_INC}:${NETCDF_INC}
	ICDF = 
else
	ICDF = -I$(MPI_INC) -I$(NETCDF_INC)
endif

	LCDF = $(NETCDF_LIB)



# Fortran loader name and options: 
 	LDR = $(FCOMP)  $(CFTFLAGS)

# Fortran compiler options/flags: (Legacy to have both CFTFLAGS & FFLAGS)
# - Optimized: (-ip = additional interprocedural optimizations for single-file compilation. )
ifeq ($(FCOMP),mpiifx)
	FFLAGS = -O3 -ipo -fp-model fast
else
	FFLAGS = -O3 -ip -fp-model fast
endif
# - uncomment for debugging:
#	FFLAGS = -g -traceback -check all -CA -CB -CS

# --------------------------------------------------------------------

ALL = partit ncjoin ncjoin_mpi extract_data_join nc3to4z zslice sample
#r2r_bry r2r_match_topo r2r_init

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

simple_xy_par_wr: simple_xy_par_wr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o simple_xy_par_wr simple_xy_par_wr.o $(LCDF)	
	
ncjoin_mpi: ncjoin_mpi.opt ncjoin_mod.o ncjoin_mpi.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin_mpi ncjoin_mod.o ncjoin_mpi.o lenstr.o $(LCDF)
	
nc3to4z: nc3to4z.o lenstr.o read_string.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o nc3to4z nc3to4z.o lenstr.o read_string.o $(LCDF)
	
mpc: mpc.F Makefile
	$(CPP) $(CPPFLAGS) -P mpc.F > ./mpc.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc ./mpc.f
	
extract_data_join: ncjoin_mod.o extract_data_join.o  lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o extract_data_join ncjoin_mod.o extract_data_join.o lenstr.o  $(LCDF)
	
particle_join: ncjoin_mod.o particle_join.o  lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o particle_join ncjoin_mod.o particle_join.o lenstr.o  $(LCDF)

R2R_INIT = r2r_init.o r2r_interp_init.o r2r_subs.o r2r_rotate.o r2r_vert_interp.o\
           set_depth.o etch_into_land.o spln2d.o tiling.o roms_read_write.o\
           read_string.o lenstr.o

r2r_init: $(R2R_INIT)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o r2r_init $(R2R_INIT) $(LCDF)

 R2R_BRY = r2r_bry.o r2r_bry_interp.o r2r_subs.o r2r_bry_rotate.o r2r_rotate.o\
           r2r_vert_interp.o set_depth.o etch_into_land.o spln2d.o tiling.o\
           def_bry_file.o roms_read_write.o read_string.o lenstr.o

r2r_bry: $(R2R_BRY)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o r2r_bry $(R2R_BRY) $(LCDF)

 ZSLICE = zslice.o set_depth.o sigma_to_z_intr.o read_add_write.o\
	  roms_read_write.o read_string.o lenstr.o tiling.o

 IPSLICE = ipslice.o set_depth.o sigma_to_z_intr.o rho_eos.o\
	   read_add_write.o roms_read_write.o read_string.o lenstr.o tiling.o

 SAMPLE = roms_read_write.o sample.o tiling.o\
          read_add_write.o read_string.o lenstr.o

R2R_MATCH_TOPO = r2r_match_topo.o r2r_interp_init.o mrg_zone_subs.o tiling.o\
                        spln2d.o roms_read_write.o read_string.o lenstr.o
                        
r2r_match_topo: $(R2R_MATCH_TOPO)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o r2r_match_topo $(R2R_MATCH_TOPO) $(LCDF)                                  
          
zslice: $(ZSLICE)
	$(LDR) $(FFLAGS) -o zslice $(ZSLICE) $(LCDF)

ipslice: $(IPSLICE)
	$(LDR) $(FFLAGS) -o ipslice $(IPSLICE) $(LCDF)

ipslice.f: zslice.F
	$(CPP) -P $(CPPFLAGS) -DIPSLICE zslice.F | ./mpc > ipslice.f

sample: $(SAMPLE)
	$(LDR) $(FFLAGS) -o sample $(SAMPLE) $(LCDF)

clean:
	/bin/rm -f *.f *.o *.a *.s *.trace *.mod

clobber: clean
	/bin/rm -f $(ALL) *.bak
	
include Make.depend
