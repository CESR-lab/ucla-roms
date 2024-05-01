# M2 Macbook (ARMAA64) version
# Universal machine independent makefile for ROMS model
#========== ======= =========== ======== === ==== =====
# Set machine dependent definitions and rules.

include Makedefs.inc

# ROMS source code files divided into groups:
#--------------------------------------------------------------------------------
# Main model driver
KRNSRC =\
	main\
	param dimensions	ocean_vars tracers\
	boundary buffer  comm_vars coupling diag_vars\
	eos_vars grid mess_buffers mixing\
	private_scratch basic_output\
	scalars scoord strings work_mod nc_read_write\
	surf_flux bulk_frc flux_frc wec_frc	tides	bgc	analytical\
	diagnostics	step2d_FB	read_inp\
	set_weights	set_scoord	init_scalars	init_arrays\
	setup_grid1	setup_grid2\
	set_nudgcof  sponge_tune roms_read_write extract_data
# 3D time stepping engine
STPSRC =\
	prsgrd	pre_step3d4S	step3d_uv1	step3d_uv2\
	step3d_t_ISO	set_depth	omega\
	visc3d_S	t3dmix_S      particles
# Boundary conditions
BCSSRC =\
	zetabc	u2dbc_im	v2dbc_im\
	u3dbc_im	v3dbc_im	w3dbc_im	t3dbc_im
# Equation of state
EOSSRC =\
	rho_eos	alfabeta
# Vertical mixing schemes
MIXSRC =\
	lmd_vmix      lmd_kpp	lmd_swr_frac
# Diagnostics and cross-checking
DIASRC =\
        diag		timers	wvlcty	grid_stiffness\
        lenstr	setup_kwds	check_srcs	check_switches1\
        get_date	ext_copy_prv2shr		check_switches2
# MPI
MPISRC=\
	mpi_setup	hidden_mpi_vars	mpi_test_arr mpi_exchanges
# Model forcing routines
FRCSRC=\
	river_frc	pipe_frc      add_git_hash\
	random_output
# Biogeochemistry
BGCSRC=\
	bgc_ecosys_vars bgc_forces bgc_param_ini	bgc_ecosys_bec2 marbl_driver\
# I/O
IOSRC = checkdims	insert_node\
	get_init	set_forces	closecdf\

  SRCS = $(KRNSRC) $(STPSRC) $(BCSSRC) $(EOSSRC) $(MIXSRC) $(DIASRC) $(MPISRC) $(FRCSRC) $(BGCSRC) $(IOSRC)
   RCS := $(addsuffix $(PPF77_ext), $(SRCS))
  OBJS := $(addsuffix $(OBJ_ext), $(SRCS))
  SBIN = roms
 LROMS = libroms.a
 LROMS2 = 2/$(LROMS)
 LROMS3 = 2/$(LROMS) 3/$(LROMS)
 LROMS4 = 2/$(LROMS) 3/$(LROMS) 4/$(LROMS)

###########
# TARGETS #
###########

####################################################################
# Machine/compiler-specific targets:
.PHONY: gnu
gnu: $(MAKE) $(SBIN)
#gnu: FFLAGS = -O3 -fallow-argument-mismatch
#	FC = mpifort
#	@if $(FC) --version | grep -qi 'GNU Fortran';then \
#         $(MAKE); \
#	else \ 
#         exit 1; \
#	fi

#.PHONY: intel
#intel: FFLAGS=


####################################################################
# ROMS targets:

# ROMS excutable (This is the first target and hence the default):
$(SBIN): $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF) $(LMPI)
	mv a.out $(SBIN)


# Multilevel libraries

$(LROMS): $(OBJS)
	/bin/rm -f $(LROMS)
	ar r $(LROMS) $(OBJS)


check_forces: check_forces$(OBJ_ext) $(SBIN)
	cp -pv check_forces$(OBJ_ext) main$(OBJ_ext) 	
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF) $(LMPI)
	mv a.out check_forces

# Everything
.PHONY: all
all: tools depend $(SBIN)

# A program to test MPI halo exchange routines.
# - ------- -- ---- --- ---- -------- ---------
MPI_TEST = mpi_test	mpi_setup	exchange	mpi_exchange4\
							mpi_exchange8WA
MPI_TEST_RCS := $(addsuffix $(PPF77_ext), $(MPI_TEST))
MPI_TEST_OBJ := $(addsuffix $(OBJ_ext)  , $(MPI_TEST))

mpi_test: $(MPI_TEST_OBJ)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpi_test $(MPI_TEST_OBJ) $(LCDF) $(LMPI)

# Auxiliary utility programs
# --------- ------- --------
 TOOLS = cppcheck srcscheck checkkwds redefs

tools: $(TOOLS)

    TMP = .

partit_t$(PPF77_ext): partit$(UPF77_ext)
	$(CPP) -P $(CPPFLAGS) -DTRANSPOSED_MPI_NODE_ORDER partit$(UPF77_ext)| mpc > partit_t$(PPF77_ext)
partit_t: partit_t$(PPF77_ext) insert_node$(OBJ_ext) lenstr$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o partit_t partit_t$(PPF77_ext) insert_node$(OBJ_ext) lenstr$(OBJ_ext) $(LCDF)

partit: partit_t partit$(OBJ_ext) insert_node$(OBJ_ext) lenstr$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o partit partit$(OBJ_ext) insert_node$(OBJ_ext) lenstr$(OBJ_ext) $(LCDF)

cppcheck: cppcheck$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck$(OBJ_ext)

srcscheck: srcscheck$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck$(OBJ_ext)

checkkwds: checkkwds$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o checkkwds checkkwds$(OBJ_ext)

redefs: redefs$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o redefs redefs$(OBJ_ext)



checkdefs: check_switches1$(UPF77_ext) setup_kwds$(UPF77_ext)

check_switches1$(UPF77_ext): cppcheck cppdefs.opt
	./cppcheck cppdefs.opt
check_srcs$(UPF77_ext): srcscheck Makefile
	./srcscheck
setup_kwds$(UPF77_ext): checkkwds read_inp$(UPF77_ext)
	./checkkwds

.PHONY: depend
depend:  # Note: cross_matrix no longer used
	@${ROMS_ROOT}/Tools-Roms/makedep.py  Make.depend
	@echo 'Updated Make.depend (dependency list)'

# Target to create tar file.
# ------ -- ------ --- -----
tarfile: clean
	tar cvf roms.tar Make* *$(UPF77_ext) *.h etalon_data.* README.*  *.in* *.mod

# Cleaning targets
# -------- -------
.PHONY: clean
clean:
	/bin/rm -rf *$(PPF77_ext) *$(OBJ_ext) *.a *.s *.mod *.trace *~ $(COMP_FILES)

.PHONY: allclean
allclean: clean
	cd 2; make -f ./Makefile clean; cd ..
	cd 3; make -f ./Makefile clean; cd ..
	cd 4; make -f ./Makefile clean; cd ..

.PHONY: clobber
clobber: clean
	/bin/rm -f check_switches1$(UPF77_ext) setup_kwds$(UPF77_ext) check_srcs$(UPF77_ext)
	/bin/rm -f $(SBIN) $(TOOLS) nsub
	/bin/rm -f core core.*
	/bin/rm -f *_genmod.mod *_genmod$(PPF90_ext)

# Automatically generated dependency list:

include Make.depend

