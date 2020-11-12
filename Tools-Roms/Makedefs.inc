# Definitions for extended memory EMT64 Intel FORTRAN compiler
# with settings for AMD Opteron processor
# =========== === ===== === ======= ======== == =====

.SUFFIXES: .o .f .F

.F.o:
	$(CPP) -P $(CPPFLAGS) $*.F | ../mpc > $*.f
	$(CFT) -c $(FFLAGS) $*.f -o $*.o
#/bin/rm -f $*.f
.F.f:
	$(CPP) -P $(CPPFLAGS) $*.F | ../mpc > $*.f
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o
.F:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<
.f:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<
.o:
	$(CFT) -o $@ $(FFLAGS) $(LDFLAGS) $<



# C-preprocessor from Intel compiler package. Note: in this case
# __IFC is a pre-defined macro, so there is no need to include it
# into CPPFLAGS.


      CPP = /lib/cpp -traditional -D__IFC
 CPPFLAGS = -I${MPI_ROOT}/include


# Compiler settings:  -fpp2 is required only if -openmp is present.
# Not having -fpp2 here just causes compiler warning (-fpp is set to
# level 2 by -openmp), but other than that has no effect.

# Switch -pc80 increases precision of floating point operation to
# 64 bits (vs. 53 bits double precision default).

# -qp  compiles and links for function profiling with  gprof(1);
#      this is the same as specifying -p or -pg.

# Setting  FFLAGS = -O2 -mp (or lower optimization level) is needed
# to pass ETALON_CHECK: -O3 causes roundoff-level differences from
# the length of innermost i-loop (the results still pass ETALON_CHECK
# if  NP_XI = NSUB_X = 1, regardless of partition in ETA-direction). 

# Flags collected under LARGE_MEM_FLAG are needed only if exceeding
# 2 GByte memory limit: both -mcmodel (sets size of addressing
# pointers) and -i-dynamic (ensures linking with proper version of
# Intel dynamic libraries) must also be specified.


# OMP_FLAG = -fpp2 -openmp

# LARGE_MEM_FLAG = -mcmodel=medium -i-dynamic

 CFTFLAGS = -pc80 -auto -stack_temps

 CFT = ifort $(OMP_FLAG) $(CFTFLAGS) $(LARGE_MEM_FLAG)

 LDR = mpif90 $(OMP_FLAG) $(CFTFLAGS) $(LARGE_MEM_FLAG)

 FFLAGS =  -O3 -IPF_fma -ip $(LARGE_MEM_FLAG)
# Debug options:
# FFLAGS = -g -check all -CA -CB -CS
# FFLAGS = -warn all -g -check arg_temp_created,bounds,pointers,uninit,format,output_conversion

 COMP_FILES = work.pc work.pcl ifort*

 LCDF = -lnetcdff

