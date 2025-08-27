.SUFFIXES: .o .f .F .f90

.F.o:
	@$(CPP) -P $*.F | mpc > $*.f
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f
	@/bin/rm -f $*.f
.F.f:
	$(CPP) -P  $*.F | mpc > $*.f
.f.o:
	$(FCOMP) -c $(FFLAGS) -o $*.o $*.f $(LCDF)

     CPP = cpp -traditional

# Path names (non-hydrostatic & NetCDF libraries):

	NHMG_ROOT = $(ROMS_ROOT)/NHMG
	NHMG_INC  = $(NHMG_ROOT)/include
	NHMG_LIB  = -L$(NHMG_ROOT)/lib -lnhmg

	ROMS_SRC  = $(ROMS_ROOT)/src

	NETCDF_INC = $(NETCDF_ROOT)/include
	NETCDF_LIB = -L$(NETCDF_ROOT)/lib -lnetcdff
	NETCDFC_LIB = -L$(NETCDFC_ROOT)/lib -lnetcdf

	MPI_INC = ${MPIHOME}/include/mpi

	LCDF = $(NETCDF_LIB) $(NETCDFC_LIB)

CPATH = ${MPI_INC}:${NETCDF_INC}

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

TSRC = test.F test_routines.F nc_read_write.F

  TCS = $(TSRC:.F=.f)
  TBJS = $(TSRC:.F=.o) 

# Everything
all: depend test

test:  $(TBJS)
	@echo 'Making: $(TBJS)'
	$(LDR) $(FFLAGS) $(LDFLAGS) -o test $(TBJS) $(LCDF)

depend:  # Note: cross_matrix no longer used
	@makedepf90 $(TSRC) > Make.depend
	@echo 'Updated Make.depend (dependency list)'

# Cleaning targets
# -------- -------
clean:
	/bin/rm -rf *.f *.o *.a *.s *.mod *.trace *~ test

# Automatically generated dependency list:

include Make.depend

