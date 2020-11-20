# Platform-specific definitions for Intel(R) Ifort(TM)


.SUFFIXES: .o .f .F .f90

.F.o:
	@$(CPP) -P $(CPPFLAGS) $*.F | ./mpc > $*.f
	$(CFT) -c $(FFLAGS) -o $*.o $*.f
	@/bin/rm -f $*.f
.F.f:
	$(CPP) -P $(CPPFLAGS) $*.F | ./mpc > $*.f
.f.o:
	$(CFT) -c $(FFLAGS) -o $*.o $*.f

.f90.o:
	$(CFT) -c $(FFLAGS) -o $*.o $*.f90
.F:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<
.f:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<
.f90:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<
.o:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<

# C-preprocessor from Intel compiler package "fpp" treats __IFC as
# a pre-defined, so there is no need to include it into CPPFLAGS, but
# the standard CPP is not aware of this.

      CPP = /lib/cpp -traditional -D__IFC
 CPPFLAGS = -I${MPI_ROOT}/include -I${NETCDF_INC}

# Compiler settings:  -fpp2 is required only if -openmp is present.
# Not having -fpp2 here just causes compiler warning (-fpp is set to
# level 2 by -openmp), but other than that has no effect.

# Switch -pc80 increases precision of floating point operation to
# 64 bits (vs. 53 bits double precision default).
#
# -qp  compiles and links for function profiling with  gprof(1);
#      this is the same as specifying -p or -pg.
#
# Setting  FFLAGS = -O2 -mp (or lower optimization level) is needed
# to pass ETALON_CHECK: -O3 causes roundoff-level differences from
# the length of innermost i-loop (the results still pass ETALON_CHECK
# if  NP_XI = NSUB_X = 1, regardless of partition in ETA-direction). 
# As of ifort v. 11.0.xxx -mp is superseeded by -fp-model <arg> flag.

# Flags collected under LARGE_MEM_FLAG are needed only if exceeding 
# 2 GBytes of memory: both -mcmodel (sets size of addressing pointers)
# and -i-dynamic (ensures proper linking with Intel dynamic libraries
# must be specified.

VPATH = ${MPI_INC}:${NETCDF_INC}:${ROMS_ROOT}/NHMG/include
CPATH = ${MPI_INC}:${NETCDF_INC}:${ROMS_ROOT}/NHMG/include

# OMP_FLAG = -fpp -openmp

# LARGE_MEM_FLAG = -mcmodel=medium -i-dynamic


#  CFTFLAGS = -pc80 -auto -stack_temps -warn unused
#  CFTFLAGS = -pc64 -auto -stack_temps -warn all
#  CFTFLAGS = -pc64 -auto -stack_temps -check all
   CFTFLAGS = -pc64 -auto -stack_temps 


     CFT = ifort $(OMP_FLAG) $(CFTFLAGS) $(LARGE_MEM_FLAG)

     LDR = mpif90 $(OMP_FLAG) $(CFTFLAGS) $(LARGE_MEM_FLAG)
 
  NHMG_ROOT = $(ROMS_ROOT)/NHMG
  NHMG_LIB  = $(NHMG_ROOT)/lib
  NHMG_INC  = $(NHMG_ROOT)/include

#  FFLAGS = -g -traceback -check all -CA -CB -CS -I./NHMG/include -I/opt/netcdf-4.6.1/include -I$(ROMS_ROOT)/ncio-master/include
#  FFLAGS = -O3 -IPF_fma -ip -I$(NHMG_INC) -I$(NETCDF_INC)
  FFLAGS = -O3 -IPF_fma -ip


  LCDF = -lnetcdff -L$(NHMG_LIB) -lnhmg
