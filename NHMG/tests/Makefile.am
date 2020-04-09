########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2016
########################################################################
SRCDIR = ../src
LIBNHMGDIR = ${SRCDIR}

AM_FCFLAGS  = ${AM_CPPFLAGS} $(FCFLAGS_f90) -I$(top_srcdir)/src
AM_FFLAGS   = ${AM_FCFLAGS}
AM_LDFLAGS  = $(TESTLIBS)
AM_LDFLAGS += ${top_builddir}/src/libnhmg.la -lnhmg

##
SUBDIRS	       = .

bin_SCRIPTS = nhmg_namelist ts_namelist rt_namelist fh_namelist cuc_namelist
EXTRA_DIST  = $(bin_SCRIPTS)

############################################################################
# Compile test only if --enable-netcdf is present in the configure options #
############################################################################
if NETCDF_CMP

if NHMG_TESTS

F90_TESTS = testseamount testseamount_rnd testrndtopo testfillhalo testcuc

check_PROGRAMS = $(F90_TESTS)
TESTS          = $(F90_TESTS)

bin_PROGRAMS   =     \
    testseamount     \
	testseamount_rnd \
	testrndtopo      \
	testfillhalo     \
	testcuc          \
    nhmg_ncjoin

##
testseamount_SOURCES=       \
	mg_zr_hz.f90            \
	mg_mpi_exchange_ijk.f90 \
    mg_setup_tests.f90      \
	mg_testseamount.f90
##
testseamount_rnd_SOURCES=   \
	mg_zr_hz.f90            \
	mg_mpi_exchange_ijk.f90 \
    mg_setup_tests.f90      \
	mg_testseamount_rnd.f90

##
testrndtopo_SOURCES=        \
	mg_zr_hz.f90            \
	mg_mpi_exchange_ijk.f90 \
    mg_setup_tests.f90      \
	mg_testrndtopo.f90

##
testfillhalo_SOURCES=       \
	mg_testfillhalo.f90

##
testcuc_SOURCES =           \
	mg_zr_hz.f90            \
	mg_mpi_exchange_ijk.f90 \
	mg_setup_tests.f90      \
	mg_testcuc.f90 

##
nhmg_ncjoin_SOURCES =       \
	nhmg_ncjoin.f90

endif

endif   

clean-local:
	-rm -f *.mod *.o

cleanf:
	-rm -f *.nc fort.* *.dump *.xml gmon.out *.log *__genmod*

##
