#include "cppdefs.h"

      subroutine put_global_atts (ncid, ierr)

      use param
#if defined SFLX_CORR && defined SALINITY
      use surf_flux, only: dSSSdt
#endif
      use eos_vars
      use ncvars
      use scalars

      implicit none
      integer ncid, ierr, nf_ftype, lvar,lenstr
#include "scoord.h"
#include "netcdf.inc"
#include "strings.h"

#if defined MPI && defined PARALLEL_FILES
      integer*4 ibuff(4), npart

! Put global attribute 'partition' which identifies subdomain
! within the processor grid individually for each file.

      npart=4
      ibuff(1)=mynode
      ibuff(2)=nnodes
      if (WESTERN_MPI_EDGE) then
        ibuff(3)=iSW_corn+iwest
      else
        ibuff(3)=iSW_corn+iwest+1
      endif
      if (SOUTHERN_MPI_EDGE) then
        ibuff(4)=jSW_corn+jsouth
      else
        ibuff(4)=jSW_corn+jsouth+1
      endif

      ierr=nf_put_att_int (ncid, nf_global, 'partition',
     &                             nf_int, npart, ibuff)
#endif

! Decide about output file type and precision for floating point
! variables (restart file always has the same precision as in the
! code, typically double precision, while all others can be made
! single precision.

      if (ncid == ncrst) then
        nf_ftype=NF_FTYPE
      else
        nf_ftype=NF_FOUT
      endif


! Title of the run

      lvar=lenstr(title)
      ierr=nf_put_att_text(ncid, nf_global, 'title',  lvar,
     &                                       title(1:lvar))

! Type of the output file

      if (ncid == ncrst) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  17,
     &                                 'ROMS restart file')
      elseif (ncid == nchis) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  17,
     &                                 'ROMS history file')
#ifdef AVERAGES
      elseif (ncid == ncavg) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  18,
     &                                'ROMS averages file')
#endif
      endif

! Grid and initial conditions file names

#ifndef ANA_GRID
      lvar=lenstr(grdname)
      ierr=nf_put_att_text(ncid, nf_global, 'grid_file',lvar,
     &                                       grdname(1:lvar))
#endif
#ifndef ANA_INITIAL
      lvar=lenstr(ininame)
      ierr=nf_put_att_text(ncid, nf_global, 'init_file',lvar,
     &                                       ininame(1:lvar))
#endif

! Time stepping parameters

      ierr=nf_put_att_int(ncid, nf_global, 'ntimes',  nf_int, 1,ntimes)
      ierr=nf_put_att_int(ncid, nf_global, 'ndtfast', nf_int,1,ndtfast)
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'dt',    nf_ftype, 1,  dt)
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'dtfast',nf_ftype,  1,
     &                                                          dtfast)
      ierr=nf_put_att_int (ncid, nf_global, 'nwrt',  nf_int, 1, nwrt)
#ifdef AVERAGES
      ierr=nf_put_att_int (ncid, nf_global, 'ntsavg',nf_int, 1, ntsavg)

      ierr=nf_put_att_int (ncid, nf_global, 'navg',  nf_int, 1,   navg)
#endif

! Vertical coordinate type, control parameters "theta_s", "theta_b",
! "hc", and stretching curves "Cs_w", "Cs_r" at W- and RHO-points.

#ifdef SOLVE3D
# ifdef VERT_COORD_TYPE_SM09
      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType', 4, 'SM09')
# elif defined VERT_COORD_TYPE_KAU
      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType', 3, 'KAU')
# elif defined VERT_COORD_TYPE_SH94
      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType', 4, 'SH94')
# endif

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'theta_s', nf_ftype, 1,
     &                                                        theta_s)
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'theta_b', nf_ftype, 1,
     &                                                        theta_b)
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'hc', nf_ftype, 1, hc)

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Cs_w',nf_ftype,N+1,Cs_w)

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Cs_r',nf_ftype, N,Cs_r)

! Boussinesq reference density and EOS coefficients for linear EOS.

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'rho0', nf_ftype, 1, rho0)
      ierr=nf_put_att_text (ncid, nf_global, 'rho0_units', 6, 'kg/m^3')
# ifndef NONLIN_EOS
      ierr=nf_put_att_FTYPE(ncid, nf_global,'Tcoef',nf_ftype, 1, Tcoef)
      ierr=nf_put_att_text (ncid, nf_global,'Tcoef_units', 13,
     &                                                 'deg Celsius-1')

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'T0',  nf_ftype, 1, T0)
      ierr=nf_put_att_text (ncid, nf_global, 'T0_units', 11,
     &                                                 'deg Celsius')

#  ifdef SALINITY
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Scoef',nf_ftype, 1,Scoef)
      ierr=nf_put_att_text (ncid, nf_global, 'Scoef_units', 5, 'PSU-1')

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'S0',  nf_ftype, 1, S0)
      ierr=nf_put_att_text (ncid, nf_global, 'S0_units',     3, 'PSU')
#  endif
# endif
#endif /* SOLVE3D */

! Horizontal viscosity and mixing coefficients.

#ifdef UV_VIS2
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'visc2',nf_ftype, 1,visc2)
      ierr=nf_put_att_text (ncid,nf_global, 'visc2_units', 5, 'm^2/s')

! Lateral slipperiness parameter

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'gamma2', nf_ftype, 1,
     &                                                       gamma2)
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'tnu2', nf_ftype, 1,tnu2)
      ierr=nf_put_att_text (ncid,nf_global, 'tnu2_units',  5,'m^2/s')
# endif
# ifdef MY25_MIXING
#  ifdef Q2_DIF2
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'q2nu2',nf_ftype, 1,q2nu2)
      ierr=nf_put_att_text (ncid,nf_global, 'q2nu2_units', 5, 'm^2/s')
#  endif
#  ifdef Q2_DIF4
      ierr=nf_put_att_FTYPE(ncid,nf_global,'q2nu4',nf_ftype, 1,q2nu4)
      ierr=nf_put_att_text (ncid,nf_global,'q2nu4_units', 5, 'm^4/s')
#  endif
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING

! Background vertical viscosity and tracer mixing coefficients.

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Akv_bak', nf_ftype, 1,
     &                                                      Akv_bak)
      ierr=nf_put_att_text (ncid, nf_global, 'Akv_bak_units', 5,
     &                                                      'm^2/s')

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Akt_bak',nf_ftype, NT,
     &                                                      Akt_bak)
      ierr=nf_put_att_text (ncid, nf_global, 'Akt_bak_units', 5,
     &                                                      'm^2/s')
# endif
# ifdef MY25_MIXING
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'Akq_bak',nf_ftype, 1,
     &                                                     Akq_bak)
      ierr=nf_put_att_text (ncid, nf_global, 'Akq_bak_units', 5,
     &                                                     'm^2/s')
# endif
#endif

! Binding velocity for external data forcing at open boundaries

#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING \
  || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M2NUDGING \
                       || defined M3NUDGING  || defined WKB_FRC_BRY

      ierr=nf_put_att_FTYPE(ncid, nf_global, 'ubind',  nf_ftype,  1,
     &                                                        ubind)
      ierr=nf_put_att_text (ncid, nf_global, 'ubind_units', 3,'m/s')
#endif

! Maximum viscosity in sponge layer

#ifdef SPONGE
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'v_sponge', nf_ftype, 1,
     &                                                     v_sponge)
      ierr=nf_put_att_text (ncid, nf_global, 'v_sponge_units',  5,
     &                                                      'm^2/s')

#endif

! Sea surface salinity correction "piston velocity"

#if defined SFLX_CORR && defined SALINITY
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'dSSSdt',  nf_ftype, 1,
     &                                        dSSSdt*(100.*day2sec))
      ierr=nf_put_att_text (ncid, nf_global, 'dSSSdt_units',  6,
     &                                                     'cm/day')
#endif

! Bottom drag parameters.

      ierr=nf_put_att_FTYPE(ncid,nf_global, 'rdrg', nf_ftype,1, rdrg)
      ierr=nf_put_att_text (ncid,nf_global, 'rdrg_units',    3,'m/s')

      ierr=nf_put_att_FTYPE(ncid,nf_global, 'rdrg2',nf_ftype,1,rdrg2)
      ierr=nf_put_att_text (ncid,nf_global, 'rdrg2_units', 14,
     &                                              'nondimensional')
#ifdef SOLVE3D
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'Zob',  nf_ftype, 1, Zob)
      ierr=nf_put_att_text (ncid,nf_global, 'Zob_units',      1, 'm')
#endif

#ifdef WKB_WWAVE

! breaking parameters for Thornton and Guza (1986),
! Church and Thornton (1993)

      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_amp',nf_ftype,1,wkb_amp)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_prd',nf_ftype,1,wkb_prd)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_ang',nf_ftype,1,wkb_ang)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_btg',nf_ftype,1,wkb_btg)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_gam',nf_ftype,1,wkb_gam)
# ifdef SURFACE_ROLLER
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_rsb',nf_ftype,1,wkb_rsb)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'wkb_roller',
     &                                            nf_ftype,1,wkb_roller)
# endif
#endif



! List of source-code file names and activated CPP-switches

      lvar=lenstr(srcs)
      ierr=nf_put_att_text (ncid, nf_global, 'SRCS', lvar,srcs(1:lvar))
      lvar=lenstr(cpps)
      ierr=nf_put_att_text (ncid, nf_global, 'CPPS', lvar,cpps(1:lvar))
      end
