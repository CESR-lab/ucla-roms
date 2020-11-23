# Universal machine independent makefile for ROMS model
#========== ======= =========== ======== === ==== =====
# Set machine dependent definitions and rules.

include Makedefs.inc
#include Makedefs.Intel_MPI
#include Makedefs.GCC

# Configuration for ROMS model:  All source code files are sorted into
# ------------- --- ---- ------  eight groups separated by blanc lines:
#   1) main driving part;          5) vertical mixing schemes;
#   2) 3D time stepping engine;    6) diagnostic and cross-checking;
#   3) boundary conditions;        7) netCDF I/O routines;
#   4) sea-water EOS routines;     8) model forcing routines;

KRNSRC = 		dimensions.F	read_write.F tracers.F\
	surf_flux.F bulk_frc.F flux_frc.F wec_frc.F	tides.F\
 step2d_FB.F	read_inp.F\
	set_weights.F	set_scoord.F	init_scalars.F	init_arrays.F\
	ana_init.F	setup_grid1.F	setup_grid2.F	visc_rescale.F\
	set_nudgcof.F	ana_grid.F	analytical.F\
\
	prsgrd.F	pre_step3d4S.F	step3d_uv1.F	step3d_uv2.F\
	step3d_t_ISO.F	set_depth.F	omega.F\
	visc3d_GP.F	t3dmix_GP.F	ana_vmix.F\
\
	zetabc.F	u2dbc_im.F	v2dbc_im.F\
	u3dbc_im.F	v3dbc_im.F	w3dbc_im.F	t3dbc_im.F\
\
	rho_eos.F	alfabeta.F\
\
	lmd_vmix.F      lmd_kpp.F	lmd_swr_frac.F	bvf_mix.F\
\
        diag.F		timers.F	wvlcty.F	grid_stiffness.F\
        lenstr.F	setup_kwds.F	check_srcs.F	check_switches1.F\
        get_date.F	ext_copy_prv2shr.F		check_switches2.F\
\
	mpi_setup.F	exchange.F	mpi_exchange8TA.F\
\
	river_frc.F	main.F		diagnostics.F

#ab_ratio.F	alfabeta.F\
#mpi_exchange4.F	mpi_exchange8TA.F\


IOSRC = checkdims.F	find_record.F	insert_node.F	ncdf_read_write.F\
	get_grid.F	get_init.F	wrt_grid.F	put_global_atts.F\
	def_rst.F	wrt_rst.F	def_his.F	wrt_his.F\
	set_avg.F	wrt_avg.F\
	get_forces.F	get_bry_all.F	set_forces.F\
	nf_read_bry.F\
	get_ssh.F	get_tclima.F	get_uclima.F	closecdf.F


#FLTSRC = init_floats.F	init_arrays_floats.F	get_init_flts.F\
#	interp_rho.F	rhs_floats.F	step_floats.F\
#	def_floats.F	wrt_floats.F

# old/obsolete/non-parallelized stuff
#       pp_vmix.F	get_wwave.F	smol_adv.F	get_bsedim.F

  SRCS = $(KRNSRC) $(IOSRC)
   RCS = $(SRCS:.F=.f)
  OBJS = $(RCS:.f=.o) 
  SBIN = roms
 LROMS = libroms.a
 LROMS2 = 2/$(LROMS)
 LROMS3 = 2/$(LROMS) 3/$(LROMS)
 LROMS4 = 2/$(LROMS) 3/$(LROMS) 4/$(LROMS)

# ROMS excutable (This is the first target and hence the default):
$(SBIN): mpc $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF) $(LMPI)
	mv a.out $(SBIN)


# Multilevel libraries

$(LROMS): $(OBJS)
	/bin/rm -f $(LROMS)
	ar r $(LROMS) $(OBJS)


check_forces: check_forces.o $(SBIN)
	cp -pv check_forces.o main.o 	
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF) $(LMPI)
	mv a.out check_forces

# Everything
all: mpc tools depend $(SBIN)

# A program to test MPI halo exchange routines.
# - ------- -- ---- --- ---- -------- ---------
MPI_TEST = mpi_test.F	mpi_setup.F	exchange.F	mpi_exchange4.F\
							mpi_exchange8WA.F
MPI_TEST_RCS = $(MPI_TEST:.F=.f)
MPI_TEST_OBJ = $(MPI_TEST_RCS:.f=.o)

mpi_test: $(MPI_TEST_OBJ)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpi_test $(MPI_TEST_OBJ) $(LCDF) $(LMPI)

# Auxiliary utility programs
# --------- ------- --------
 TOOLS = mpc cppcheck srcscheck checkkwds redefs

tools: $(TOOLS)

    TMP = .

mpc: mpc.F Makefile
	$(CPP) $(CPPFLAGS) -P mpc.F > $(TMP)/mpc.f
	$(CFT) $(FFLAGS) $(LDFLAGS) -o mpc $(TMP)/mpc.f
#	/bin/rm -rf $(TMP)/mpc.f

partit_t.f: partit.F
	$(CPP) -P $(CPPFLAGS) -DTRANSPOSED_MPI_NODE_ORDER partit.F| mpc > partit_t.f
partit_t: partit_t.f insert_node.o lenstr.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o partit_t partit_t.f insert_node.o lenstr.o $(LCDF)

partit: partit_t partit.o insert_node.o lenstr.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o lenstr.o $(LCDF)

cppcheck: cppcheck.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck.o

srcscheck: srcscheck.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck.o

checkkwds: checkkwds.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o checkkwds checkkwds.o

redefs: redefs.o
	$(CFT) $(FFLAGS) $(LDFLAGS) -o redefs redefs.o



checkdefs: check_switches1.F setup_kwds.F

check_switches1.F: cppcheck cppdefs.h
	./cppcheck cppdefs.h
check_srcs.F: srcscheck Makefile
	./srcscheck
setup_kwds.F: checkkwds read_inp.F
	./checkkwds

depend:  # Note: cross_matrix no longer used
	@${ROMS_ROOT}/Tools-Roms/makedepf90 $(SRCS) > Make.depend
	@echo 'Updated Make.depend (dependency list)'

plotter: plotter.o
	$(CFT) -o plotter plotter.o $(LIBNCAR)

# Target to create tar file.
# ------ -- ------ --- -----
tarfile: clean
	tar cvf roms.tar Make* *.F *.h etalon_data.* README.*  *.in* *.mod

# Cleaning targets
# -------- -------
clean:
	/bin/rm -rf *.f *.o *.a *.s *.mod *.trace *~ $(COMP_FILES)
allclean: clean
	cd 2; make -f ./Makefile clean; cd ..
	cd 3; make -f ./Makefile clean; cd ..
	cd 4; make -f ./Makefile clean; cd ..

clobber: clean
	/bin/rm -f check_switches1.F setup_kwds.F check_srcs.F
	/bin/rm -f $(SBIN) $(TOOLS) nsub
	/bin/rm -f core core.*
	/bin/rm -f *_genmod.mod *_genmod.f90

# Automatically generated dependency list:

include Make.depend
