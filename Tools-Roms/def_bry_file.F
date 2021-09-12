      subroutine def_bry_file(fname, ncid, xi_rho,eta_rho,N,
     &                OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH,
     &                        theta_s, theta_b, hc, Cs_w,Cs_r)
      implicit none
      character(len=*) fname
      integer  ncid, xi_rho,eta_rho,N
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH
      real(kind=8) theta_s,theta_b,hc, Cs_w(0:N),Cs_r(N)
      integer r2d_EW(2), u2d_EW(2), r3d_EW(3), u3d_EW(3),
     &        r2d_NS(2), v2d_NS(2), r3d_NS(3), v3d_NS(3),
     &        vert_dim, time_dim, varid, ierr,
     &        old_fill_mode, lenstr
#include "spval.h"
#include "netcdf.inc"

      write(*,'(1x,3A)',advance='no') 'def_bry_file :: creating ''',
     &                            fname(1:lenstr(fname)),  ''' ...'

      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)

! Dimensions

      ierr=nf_def_dim(ncid, 'xi_rho',  xi_rho,      r2d_EW(1))
      ierr=nf_def_dim(ncid, 'xi_u',    xi_rho-1,    u2d_EW(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho,     r2d_NS(1))
      ierr=nf_def_dim(ncid, 'eta_v',   eta_rho-1,   v2d_NS(1))
      ierr=nf_def_dim(ncid, 's_rho',   N,            vert_dim)
      ierr=nf_def_dim(ncid, 'rec_time',nf_unlimited, time_dim)

      r3d_EW(1)=r2d_EW(1) ; r2d_EW(2)=time_dim
      u3d_EW(1)=u2d_EW(1) ; u2d_EW(2)=time_dim
      r3d_NS(1)=r2d_NS(1) ; r2d_NS(2)=time_dim
      v3d_NS(1)=v2d_NS(1) ; v2d_NS(2)=time_dim

      r3d_EW(2)=vert_dim  ; r3d_EW(3)=time_dim
      u3d_EW(2)=vert_dim  ; u3d_EW(3)=time_dim
      r3d_NS(2)=vert_dim  ; r3d_NS(3)=time_dim
      v3d_NS(2)=vert_dim  ; v3d_NS(3)=time_dim


! Recoord S-coordinate control parameters "theta_s", "theta_b", "hc",
! and stretching curves "Cs_w", "Cs_r" at vertical W- and RHO-points.

      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType',4,'SM09')

      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double, 1,
     &                                                       theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double, 1,
     &                                                       theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc', nf_double, 1, hc)

      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double, N+1,
     &                                                          Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double, N,
     &                                                         Cs_r)

! Time

      ierr=nf_def_var (ncid, 'bry_time', nf_double, 1, time_dim, varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')

! Side boundary forcing variables

      if (OBC_WEST) then
        ierr=nf_def_var(ncid, 'zeta_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'ubar_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,  1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'vbar_west',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'u_west',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'v_west',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'temp_west',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)

        ierr=nf_def_var(ncid, 'salt_west',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif

      if (OBC_EAST) then
        ierr=nf_def_var(ncid, 'zeta_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'ubar_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'vbar_east',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'u_east',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'v_east',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)

        ierr=nf_def_var(ncid, 'temp_east',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)

        ierr=nf_def_var(ncid, 'salt_east',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif

      if (OBC_SOUTH) then
        ierr=nf_def_var(ncid, 'zeta_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'ubar_south',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)

        ierr=nf_def_var(ncid, 'vbar_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'u_south',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'v_south',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'temp_south',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'salt_south',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif

      if (OBC_NORTH) then
        ierr=nf_def_var(ncid, 'zeta_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'ubar_north',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'vbar_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'u_north',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)

        ierr=nf_def_var(ncid, 'v_north',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)

        ierr=nf_def_var(ncid, 'temp_north',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)

        ierr=nf_def_var(ncid, 'salt_north',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif
      write(*,*) '...done.'
      end

