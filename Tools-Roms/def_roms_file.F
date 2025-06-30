      subroutine def_roms_file(ncid, fname, xi_rho,eta_rho,s_rho,
     &                            theta_s,theta_b, hc, Cs_w,Cs_r)

! Set up netCDF structure for a file suitable to be ROMS initial or
! climatological input file.  Note that the file itself is expected to
! be created externally - incoming argument "ncid" is a valid netCDF ID
! of a writable file, while this routine only takes care about creating
! dimensions, variables, and writing attributes.  Similarly, it does
! not call nf_enddef to finish the definitions and switch into input
! mode leaving the calling program some room to customize the
! definitions without resorting to nf_redef (which may require moving
! the data).

      implicit none
      character(len=*) fname
      integer ncid, xi_rho,xi_u, eta_rho,eta_v, s_rho, ierr,
     &                        old_fill_mode, varid,  lfnm, lenstr
      real(kind=8) theta_s,theta_b, hc, Cs_w(0:s_rho),Cs_r(s_rho)
#ifdef WITH_REC_DIM
      integer timedim
      integer, parameter :: n2d=3, n3d=4
#else
      integer, parameter :: n2d=2, n3d=3
#endif
      integer r2dgrd(n2d), u2dgrd(n2d), v2dgrd(n2d),
     &        r3dgrd(n3d), u3dgrd(n3d), v3dgrd(n3d)
      real(kind=8), parameter :: cycle_length=360.D0
#include "spval.h"
#include "netcdf.inc"

      lfnm=lenstr(fname) ;  write(*,'(1x,3A)', advance='no')
     &    'creating netCDF file ''', fname(1:lfnm), ''' ...'

      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)
      if (ierr /= nf_noerr) then
        write(*,'(/1x,A/)') '### WARNING: Cannot set nofill mode.'
      endif

! Define dimensions

      xi_u=xi_rho-1 ; eta_v=eta_rho-1
      ierr=nf_def_dim(ncid, 'xi_rho',  xi_rho,    r2dgrd(1))
      ierr=nf_def_dim(ncid, 'xi_u',    xi_u,      u2dgrd(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho,   r2dgrd(2))
      ierr=nf_def_dim(ncid, 'eta_v',   eta_v,     v2dgrd(2))
      ierr=nf_def_dim(ncid, 's_rho',   s_rho,     r3dgrd(3))
#ifdef WITH_REC_DIM
      ierr=nf_def_dim(ncid, 'time',   nf_unlimited, timedim)
#endif

      v2dgrd(1)=r2dgrd(1) ; u2dgrd(2)=r2dgrd(2)

      r3dgrd(1)=r2dgrd(1) ; r3dgrd(2)=r2dgrd(2)
      u3dgrd(1)=u2dgrd(1) ; u3dgrd(2)=u2dgrd(2)
      v3dgrd(1)=v2dgrd(1) ; v3dgrd(2)=v2dgrd(2)

      u3dgrd(3)=r3dgrd(3) ; v3dgrd(3)=r3dgrd(3)

#ifdef WITH_REC_DIM
      r2dgrd(3)=timedim ; u2dgrd(3)=timedim ; v2dgrd(3)=timedim
      r3dgrd(4)=timedim ; u3dgrd(4)=timedim ; v3dgrd(4)=timedim
#endif

! Recoord S-coordinate control parameters "theta_s", "theta_b", "hc",
! and stretching curves "Cs_w", "Cs_r" at vertical W- and RHO-points.

      ierr=nf_put_att_text(ncid,nf_global,'VertCoordType',4,'SM09')

      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double,
     &                                                 1, theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double,
     &                                                 1, theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc',nf_double,1,hc)

      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double,
     &                                              s_rho+1, Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double,
     &                                                s_rho, Cs_r)

! Time.

#ifdef SEPARATE_TIME_VARS
      ierr=nf_def_var (ncid, 'ssh_time',  nf_double, n2d-2,
     &                                     r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_put_att_double (ncid, varid, 'cycle_length',
     &                             nf_double, 1, cycle_length)

      ierr=nf_def_var (ncid, 'ssh_time',  nf_double, n2d-2,
     &                                       r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_put_att_double (ncid, varid, 'cycle_length',
     &                               nf_double, 1, cycle_length)

      ierr=nf_def_var (ncid, 'uclm_time',  nf_double, n2d-2,
     &                                      r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_put_att_double (ncid, varid, 'cycle_length',
     &                              nf_double, 1, cycle_length)

      ierr=nf_def_var (ncid, 'tclm_time',  nf_double, n2d-2,
     &                                       r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_put_att_double (ncid, varid, 'cycle_length',
     &                               nf_double, 1, cycle_length)

      ierr=nf_def_var (ncid, 'sclm_time',  nf_double, n2d-2,
     &                                       r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_put_att_double (ncid, varid, 'cycle_length',
     &                               nf_double, 1, cycle_length)
#else
      ierr=nf_def_var (ncid, 'ocean_time',  nf_double, n2d-2,
     &                                      r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name', 22,
     &                             'time since initialization')
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
#endif

! Free-surface.

      ierr=nf_def_var(ncid, 'zeta', nf_real, n2d, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                  'free-surface elevation')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                 FillValue)


! 2D momenta in XI- and ETA-directions.

      ierr=nf_def_var(ncid, 'ubar', nf_real, n2d, u2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                   'barotropic XI-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)


      ierr=nf_def_var(ncid, 'vbar', nf_real, n2d, v2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 23,
     &                                  'barotropic ETA-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

! 3D momenta in XI- and ETA-directions.

      ierr=nf_def_var(ncid, 'u', nf_real, n3d, u3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 21,
     &                                    'XI-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)


      ierr=nf_def_var(ncid, 'v', nf_real, n3d, v3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                   'ETA-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

! Tracer variables.

      ierr=nf_def_var(ncid, 'temp', nf_real, n3d, r3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 21,
     &                                   'potential temperature')
      ierr=nf_put_att_text(ncid, varid, 'units', 7,    'Celsius')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     &                                              1, FillValue)

      ierr=nf_def_var(ncid, 'salt', nf_real, n3d, r3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 8,'salinity')
      ierr=nf_put_att_text(ncid, varid, 'units',     3,     'PSU')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     &                                              1, FillValue)
      end
