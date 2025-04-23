# Makefile for UCLA ROMS model
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

#
.SUFFIXES: .o .f .F .f90

.F.o:
	@$(CPP) -P $*.F | mpc > $*.f
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f
	@/bin/rm -f $*.f
.F.f:
	$(CPP) -P  $*.F | mpc > $*.f
.f.o:
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f $(LCDF)
.f90.o:
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f90

     CPP = cpp -traditional

# Path names (non-hydrostatic & NetCDF libraries):

	NHMG_ROOT = $(ROMS_ROOT)/NHMG
	NHMG_INC  = $(NHMG_ROOT)/include
	NHMG_LIB  = -L$(NHMG_ROOT)/lib -lnhmg

	NETCDF_INC = $(NETCDF_ROOT)/include
	NETCDF_LIB = -L$(NETCDF_ROOT)/lib -lnetcdff
	NETCDFC_LIB = -L$(NETCDFC_ROOT)/lib -lnetcdf

	MPI_INC = ${MPIHOME}/include/mpi

	LCDF = $(NHMG_LIB) $(NETCDF_LIB) $(NETCDFC_LIB)

	CPATH = ${MPI_INC}:${NHMG_INC}:${NETCDF_INC}

# MARBL
#       MARBL_INC = -I${MARBL_ROOT}/include/intel-mpi
#   	MARBL_LIB = -lmarbl-intel-mpi -L${MARBL_ROOT}/lib


# Large memory runs (e.g. for bgc):
#	LARGE_MEM_FLAG = -mcmodel=medium
#	LARGE_MEM_FLAG = -mcmodel=large

# Fortran compiler flags:
#	CFTFLAGS = -pc64 -auto
	CFTFLAGS = 

# Fortran loader name and options: 
 	LDR = $(FCOMP)  $(CFTFLAGS)

# Fortran compiler options/flags: (Legacy to have both CFTFLAGS & FFLAGS)
# - Optimized: (-ip = additional interprocedural optimizations for single-file compilation. )
	FFLAGS = -O3 -ipo -fp-model fast
# - uncomment for debugging:
#	FFLAGS = -g -traceback -check all -CA -CB -CS
#	FFLAGS = -g -traceback -CA -CB -CS
# - uncomment for traceback only (traceback not as accurate as full debug mode...)
#	FFLAGS = -O3 -ip -traceback
# - for code_check script ensure consistency with -fp-model strict:
#	FFLAGS = -O3 -fp-model strict
# - for profiling with vtune to see source code, also comment '@/bin/rm -f $*.f':
#	FFLAGS = -O3 -g -debug inline-debug-info -parallel-source-info=2
# - for profiling with grof to see source code, also comment '@/bin/rm -f $*.f':
#	FFLAGS = -O3 -ip -pg


KRNSRC = 		param.F dimensions.F	ocean_vars.F tracers.F\
	boundary.F buffer.F  comm_vars.F coupling.F diag_vars.F\
	eos_vars.F grid.F hidden_mpi_vars.F mess_buffers.F mixing.F\
	mpi_test_arr.F private_scratch.F basic_output.F\
	scalars.F scoord.F strings.F work_mod.F nc_read_write.F\
	surf_flux.F bulk_frc.F flux_frc.F wec_frc.F	tides.F	bgc.F	analytical.F\
	diagnostics.F	step2d_FB.F	read_inp.F\
	set_weights.F	set_scoord.F	init_scalars.F	init_arrays.F\
	setup_grid1.F	setup_grid2.F\
	set_nudgcof.F  sponge_tune.F roms_read_write.F extract_data.F \
\
	prsgrd.F	pre_step3d4S.F	step3d_uv1.F	step3d_uv2.F\
	step3d_t_ISO.F	set_depth.F	omega.F\
	visc3d_S.F	t3dmix_S.F      particles.F	visc2d.F\
\
	zetabc.F	u2dbc_im.F	v2dbc_im.F\
	u3dbc_im.F	v3dbc_im.F	w3dbc_im.F	t3dbc_im.F\
\
	rho_eos.F	alfabeta.F\
\
	lmd_vmix.F      lmd_kpp.F	lmd_swr_frac.F	\
\
        diag.F		timers.F	wvlcty.F	grid_stiffness.F\
        lenstr.F	setup_kwds.F	check_srcs.F	check_switches1.F\
        get_date.F	ext_copy_prv2shr.F		check_switches2.F\
\
	mpi_setup.F\
\
	river_frc.F	pipe_frc.F      add_git_hash.F\
	random_output.F cstar_output.F mpi_exchanges.F	zslice_output.F\
\
	bgc_ecosys_vars.F bgc_forces.F bgc_param_ini.F	bgc_ecosys_bec2.F marbl_driver.F\
\
	main.F	

IOSRC = checkdims.F	insert_node.F\
	get_init.F	set_forces.F	closecdf.F\

  SRCS = $(KRNSRC) $(IOSRC)
   RCS = $(SRCS:.F=.f)
  OBJS = $(RCS:.f=.o) 
  SBIN = roms

# ROMS excutable (This is the first target and hence the default):
$(SBIN): $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF)
	mv a.out $(SBIN)

# Everything
all: tools depend $(SBIN)

# Auxiliary utility programs
# --------- ------- --------
 TOOLS = cppcheck srcscheck checkkwds redefs

tools: $(TOOLS)

    TMP = .

cppcheck: cppcheck.o
	$(FCOMP) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck.o

srcscheck: srcscheck.o
	$(FCOMP) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck.o

checkkwds: checkkwds.o
	$(FCOMP) $(FFLAGS) $(LDFLAGS) -o checkkwds checkkwds.o


checkdefs: check_switches1.F setup_kwds.F

check_switches1.F: cppcheck cppdefs.opt
	./cppcheck cppdefs.opt
check_srcs.F: srcscheck Makefile
	./srcscheck
setup_kwds.F: checkkwds read_inp.F
	./checkkwds

depend:  # Note: cross_matrix no longer used
	@${ROMS_ROOT}/Tools-Roms/makedepf90 $(SRCS) > Make.depend
	@echo 'Updated Make.depend (dependency list)'

# Target to create tar file.
# ------ -- ------ --- -----
tarfile: clean
	tar cvf roms.tar Make* *.F *.h etalon_data.* README.*  *.in* *.mod

# Cleaning targets
# -------- -------
clean:
	/bin/rm -rf *.f *.o *.a *.s *.mod *.trace *~ $(COMP_FILES)

# Automatically generated dependency list:

include Make.depend

