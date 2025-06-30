#define MASK_LAND_DATA

      subroutine set_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer nx,ny, i,j
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        do i=1,nx                        ! Convert RHO-point mask,
          if (rmsk(i,j) < 0.5D0) then    ! real(kind=8) array into
            mask(i,j)=0                  ! short 2-byte integer.
          else
            mask(i,j)=1                  ! Note that here U- and V-
          endif                          ! masks are set to zero only
        enddo                            ! for U- and V-point points
        do i=2,nx
          if (rmsk(i,j) < 0.5D0 .and. rmsk(i-1,j) < 0.5D0) then
            umask(i,j)=0
          else                           ! which are FULLY INSIDE land
            umask(i,j)=1                 ! (hence ".and." logic instead
          endif                          ! of ".or." in the actual ROMS
        enddo                            ! model), while points of
        if (j > 1) then                  ! normal velocity components
          do i=1,nx
            if (rmsk(i,j) < 0.5D0 .and. rmsk(i,j-1) < 0.5D0) then
              vmask(i,j)=0
            else                         ! on the coast line - the
              vmask(i,j)=1               ! ones where normal velocities
            endif                        ! must be set to zero due to
          enddo                          ! no-flux boundary conditions
        endif                            ! - are considered to be part
      enddo                              ! of the the solution and
      end                                ! therefore left unmasked.


      subroutine set_contour_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer nx,ny, i,j, iw2,iw,ie,js2,js,jn
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        js2=max(1,j-2) ; js=max(1,j-1) ; jn=min(j+1,ny)
        do i=1,nx
          iw=max(1,i-1) ; ie=min(i+1,nx)
          if (rmsk(i,j) < 0.5D0 .and. ( rmsk(ie,j) > 0.5D0
     &       .or. rmsk(i,jn) > 0.5D0 .or. rmsk(iw,j) > 0.5D0
     &                       .or. rmsk(i,js) > 0.5D0 )) then
            mask(i,j)=0
          else
            mask(i,j)=1
          endif
        enddo
        do i=2,nx
          iw2=max(1,i-2) ; iw=max(1,i-1) ; ie=min(i+1,nx)
          if ( rmsk(i-1,j) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     &        (  rmsk(i-1,jn) > 0.5D0 .or. rmsk(i,jn) > 0.5D0  .or.
     &       rmsk(iw2,j) > 0.5D0      .or.     rmsk(ie,j) > 0.5D0 .or.
     &           rmsk(i-1,js) > 0.5D0 .or. rmsk(i,js) > 0.5D0 )) then
            umask(i,j)=0
          else
            umask(i,j)=1
          endif
        enddo
        if (j > 1) then
          do i=1,nx
            iw=max(1,i-1) ; ie=min(i+1,nx)
            if ( rmsk(i,j-1) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     &          (             rmsk(i,jn) > 0.5D0 .or.
     &           rmsk(iw,j  ) > 0.5D0 .or. rmsk(ie,j  ) > 0.5D0 .or.
     &           rmsk(iw,j-1) > 0.5D0 .or. rmsk(ie,j-1) > 0.5D0 .or.
     &                        rmsk(i,js2) > 0.5D0 )) then
              vmask(i,j)=0
            else
              vmask(i,j)=1
            endif
          enddo
        endif
      enddo
      end

      subroutine apply_mask(nx,ny,N, mask,q)
      integer nx,ny,N, i,j,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(nx,ny,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#else
# define FillValue 0.
#endif
      if (N > 1) then
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) then
              do k=1,N
                q(i,j,k)=FillValue
              enddo
            endif
          enddo
        enddo
      else
        k=1
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) q(i,j,k)=FillValue
          enddo
        enddo
      endif
      end

      subroutine bry_apply_mask(nx,N, mask,q)
      integer nx,N, i,k
      integer(kind=2) mask(nx)    ! Apply land mask
      real (kind=4) q(nx,N)       ! to horizontally
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      if (N > 1) then             ! one-dimensional
        do i=1,nx                 ! array. Note that
          if (mask(i) == 0) then  ! mask array itself
            do k=1,N              ! is horizontally
              q(i,k)=FillValue    ! one-dimensional
            enddo                 ! too.
          endif
        enddo
      else
        k=1
        do i=1,nx
          if (mask(i) == 0) q(i,k)=FillValue
        enddo
      endif
      end




/* ------> OBSOLETE ------->
      subroutine apply_mask_west(nx,ny,N, mask,q)
      integer nx,ny,N, j,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(ny,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      if (N > 1) then
        do j=1,ny
          if (mask(1,j) == 0) then
            do k=1,N
              q(j,k)=FillValue
            enddo
          endif
        enddo
      else
        k=1
        do j=1,ny
          if (mask(1,j) == 0) q(j,k)=FillValue
        enddo
      endif
      end
      subroutine apply_mask_east(nx,ny,N, mask,q)
      integer nx,ny,N, j,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(ny,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      if (N > 1) then
        do j=1,ny
          if (mask(nx,j) == 0) then
            do k=1,N
              q(j,k)=FillValue
            enddo
          endif
        enddo
      else
        k=1
        do j=1,ny
          if (mask(nx,j) == 0) q(j,k)=FillValue
        enddo
      endif
      end
      subroutine apply_mask_south(nx,ny,N, mask,q)
      integer nx,ny,N, i,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(nx,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      if (N > 1) then
        do i=1,nx
          if (mask(i,1) == 0) then
            do k=1,N
              q(i,k)=FillValue
            enddo
          endif
        enddo
      else
        k=1
        do i=1,nx
          if (mask(i,1) == 0) q(i,k)=FillValue
        enddo
      endif
      end
      subroutine apply_mask_north(nx,ny,N, mask,q)
      integer nx,ny,N, i,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(nx,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      if (N > 1) then
        do i=1,nx
          if (mask(i,ny) == 0) then
            do k=1,N
              q(i,k)=FillValue
            enddo
          endif
        enddo
      else
        k=1
        do i=1,nx
          if (mask(i,ny) == 0) q(i,k)=FillValue
        enddo
      endif
      end
<---------------- OBSOLETE <-------  */


      subroutine compute_min_max(nx,ny, x, xmin,xmax)
      implicit none
      integer nx,ny, i,j
      real(kind=8) x(nx,ny), xmin,xmax
      xmin=x(1,1) ; xmax=x(1,1)
      do j=1,ny                         ! Compute min,max values of
        do i=1,nx                       ! array "x" irrespective of
          if (x(i,j) < xmin) then       ! what is its meaning.
            xmin=x(i,j)
          elseif (x(i,j) > xmax) then
            xmax=x(i,j)
          endif
        enddo
      enddo
      write(*,*) 'min,max =', xmin,xmax
      end

      subroutine adjust_lon_into_range(nx,ny, x, west,east)
      implicit none
      integer nx,ny, i,j
      real(kind=8) x(nx,ny), west,east, xmin,xmax

      xmin=x(1,1) ; xmax=x(1,1)         ! Compute min,max values of
      do j=1,ny                         ! longitude (array "x") whose
        do i=1,nx                       ! definition is presumed to be
          if (x(i,j) < xmin) then       ! unknown in advance, e.g., it
            xmin=x(i,j)                 ! may be from 0 to 360 or from
          elseif (x(i,j) > xmax) then   ! -180 to +180 or whatever --
            xmax=x(i,j)                 ! it is however presumed to be
          endif                         ! continuous within the array.
        enddo                           ! Check whether there is an
      enddo                             ! overlap with the externally
      write(*,*) 'min,max =', xmin,xmax ! specified west,east range.
      if (xmin > east) then             ! Adjust by adding/subtracting
        do j=1,ny                       ! 360 degrees is necessary to
          do i=1,nx                     ! make it overlap.
            x(i,j)=x(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin-360.D0,xmax-360.D0
      elseif (xmax < west) then
        do j=1,ny
          do i=1,nx
            x(i,j)=x(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin+360.D0, xmax+360.D0
      else
        write(*,*) 'no need to adjust parent-grid longitude'
      endif
      end

! The following routines find bounding indices on parent grid which
! define minimal logically rectangular patch fully containing unmasked
! portion of child grid. It is almost the same finding min,max for two
! integer arrays, expect that there are two caveats:
!
!  (i) some of indices ip,jp may be set to non-positive values to
!      indicate that that portion of child grid is outside the parent,
!      so interpolation is not possible there -- it is still OK as
!      long as these areas are under land mask -- these points are
!      ignored for the purpose of search below.  For this reason the
!      search is done by full 2D-sweep rather than moving along the
!      perimeter of child grid (which may not be even closed because
!      of having special-valued (ip,jp)-s; another good reason for
!      needing 2D-sweep is land mask; and
!
! (ii) indices [ip(i,j),jp(i,j)] mean that child-grid point i,j is
!      located somewhere within the parent-grid area bounded by 4
!      vertices
!
!        (xp(ip,jp+1),yp(ip,jp+1)) --- (xp(ip+1,jp+1),yp(ip+1,jp+1))
!                     |                              |
!                     |                              |
!            (xp(ip,jp),yp(ip,jp))  ---   (xp(ip+1,jp),yp(ip+1,jp))
!
!      so +1 is added to both ieast and jnorth at the very end.


      subroutine compute_index_bounds(ncx,ncy, ip,jp, mask, imin,imax,
     &                                                      jmin,jmax)
      implicit none
      integer ncx,ncy,ip(ncx,ncy),jp(ncx,ncy), imin,imax,jmin,jmax,i,j
      integer(kind=2) mask(ncx,ncy)
      imin=10000000 ; imax=-1  !<-- initialize to unrealistic
      jmin=10000000 ; jmax=-1  !<-- values outside the range
      do j=1,ncy
        do i=1,ncx
          if (mask(i,j) > 0 .and. ip(i,j) > 0 .and. jp(i,j) > 0) then
            imin=min(imin,ip(i,j)) ; jmin=min(jmin,jp(i,j))
            imax=max(imax,ip(i,j)) ; jmax=max(jmax,jp(i,j))
          endif
        enddo
      enddo
      imax=imax+1 ; jmax=jmax+1   !<-- because of (ii) above
      end


! Same as above, but to be applied along a 1D line. Note that this
! time imin,imax,jmin,jmax are expected to be initialized externally
! and there is no adding +1 at the end: this is because consecutive
! calls to this routine applied to the different sides of perimeter
! of the grid (e.g., open boundaries) contribute to finding the same
! bounds and it needs to be summarized externally after all done.

      subroutine r2r_bry_index_bounds(ncx, ip,jp, mask, imin,imax,
     &                                                  jmin,jmax)
      implicit none
      integer ncx, ip(ncx),jp(ncx), imin,imax,jmin,jmax, ic
      integer(kind=2) mask(ncx)
      do ic=1,ncx
        if (mask(ic) > 0 .and. ip(ic) > 0 .and. jp(ic) > 0) then
          imin=min(imin,ip(ic)) ; jmin=min(jmin,jp(ic))
          imax=max(imax,ip(ic)) ; jmax=max(jmax,jp(ic))
        endif
      enddo
      end




/*-------> OBSOLETE ------->
! The following four routine are same as above, but for 1D searches
! along one of the open boundaries.   Unlike the above, they neither
! initialize nor finalize imin,imax,jmin,jmax  leaving it to calling
! program.  They differ from each other only by having  mask(1,j) vs.
! mask(ncx,j) vs. mask(i,1) vs. mask(i,ncy), and by using index i
! instead of j as appropriate.

      subroutine comp_index_bnds_west(ncx,ncy, ip,jp, mask,
     &                                 imin,imax,jmin,jmax)
      implicit none
      integer ncx,ncy, ip(ncx),jp(ncx), imin,imax,jmin,jmax, j
      integer(kind=2) mask(ncx,ncy)
      do j=1,ncy
        if (mask(1,j) > 0 .and. ip(j) > 0 .and. jp(j) > 0) then
          imin=min(imin,ip(j)) ; jmin=min(jmin,jp(j))
          imax=max(imax,ip(j)) ; jmax=max(jmax,jp(j))
        endif
      enddo
      end
      subroutine comp_index_bnds_east(ncx,ncy, ip,jp, mask,
     &                                 imin,imax,jmin,jmax)
      implicit none
      integer ncx,ncy, ip(ncx),jp(ncx), imin,imax,jmin,jmax, j
      integer(kind=2) mask(ncx,ncy)
      do j=1,ncy
        if (mask(ncx,j) > 0 .and. ip(j) > 0 .and. jp(j) > 0) then
          imin=min(imin,ip(j)) ; jmin=min(jmin,jp(j))
          imax=max(imax,ip(j)) ; jmax=max(jmax,jp(j))
        endif
      enddo
      end

      subroutine comp_index_bnds_south(ncx,ncy, ip,jp, mask,
     &                                  imin,imax,jmin,jmax)
      implicit none
      integer ncx,ncy, ip(ncx),jp(ncx), imin,imax,jmin,jmax, i
      integer(kind=2) mask(ncx,ncy)
      do i=1,ncx
        if (mask(i,1) > 0 .and. ip(i) > 0 .and. jp(i) > 0) then
          imin=min(imin,ip(i)) ; jmin=min(jmin,jp(i))
          imax=max(imax,ip(i)) ; jmax=max(jmax,jp(i))
        endif
      enddo
      end

      subroutine comp_index_bnds_north(ncx,ncy, ip,jp, mask,
     &                                  imin,imax,jmin,jmax)
      implicit none
      integer ncx,ncy, ip(ncx),jp(ncx), imin,imax,jmin,jmax, i
      integer(kind=2) mask(ncx,ncy)
      do i=1,ncx
        if (mask(i,ncy) > 0 .and. ip(i) > 0 .and. jp(i) > 0) then
          imin=min(imin,ip(i)) ; jmin=min(jmin,jp(i))
          imax=max(imax,ip(i)) ; jmax=max(jmax,jp(i))
        endif
      enddo
      end
<---------------- OBSOLETE  <------ */



! The following routine is called from "r2r_init" and writes out all
! the incoming 2D arrays into a special-purpose netCDF file it creates
! internally. This is needed only for diagnostic purposes only and does
! not affect outcome of any computations.  Consequently the routine is
! written in a relaxed way without performing checks for netCDF errors.


      subroutine r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     &                     ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     &                              csA,snA, h, hprnt, kprnt)
      implicit none
      integer ncx,ncy, N
      integer(kind=4), dimension(ncx,ncy) :: ip,jp, ipu,jpu, ipv,jpv
      real(kind=8), dimension(ncx,ncy)  :: csA,snA, h,hprnt, xi,eta,
     &                                             xiu,etau, xiv,etav
      real(kind=8), dimension(ncx,ncy,N) :: kprnt
      integer ncid, ierr,i,j,dimids(4), ipvar, jpvar,  xivar, etavar,
     &            old_fill_mode, csvar, ipuvar,jpuvar, xiuvar,etauvar,
     &                           snvar, ipvvar,jpvvar, xivvar,etavvar,
     &                                  hvar,  hpvar,  dhvar, kpvar
      real(kind=8), allocatable, dimension(:,:) :: dh

#include "netcdf.inc"

      write(*,'(1x,2A)',advance='no') 'creating diagnostic file ',
     &                                  '''r2r_init_diag.nc''...'
      ierr=nf_create('r2r_init_diag.nc', nf_64bit_offset, ncid)
      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)

      ierr=nf_def_dim(ncid, 'xi_rho',  ncx, dimids(1))
      ierr=nf_def_dim(ncid, 'eta_rho', ncy, dimids(2))
      ierr=nf_def_dim(ncid, 's_rho',   N,   dimids(3))
!>
      ierr=nf_def_var(ncid, 'ip', nf_int, 2, dimids, ipvar)
      ierr=nf_put_att_text(ncid, ipvar, 'long_name', 37,
     &                    'parent-to-child interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jp', nf_int, 2, dimids, jpvar)
      ierr=nf_put_att_text(ncid, jpvar, 'long_name', 37,
     &                    'parent-to-child interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xi',  nf_double, 2, dimids,  xivar)
      ierr=nf_put_att_text(ncid, xivar, 'long_name', 38,
     &                  'parent-to-child interpolation offset X')
      ierr=nf_put_att_double(ncid, xivar, '_FillValue', nf_double,
     &                                                  1, -1.D0)

      ierr=nf_def_var(ncid, 'eta', nf_double, 2, dimids, etavar)
      ierr=nf_put_att_text(ncid, etavar, 'long_name', 38,
     &                  'parent-to-child interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavar,'_FillValue', nf_double,
     &                                                  1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'ipu', nf_int, 2, dimids, ipuvar)
      ierr=nf_put_att_text(ncid, ipuvar, 'long_name', 39,
     &                   'parent-to-child U-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipuvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jpu', nf_int, 2, dimids, jpuvar)
      ierr=nf_put_att_text(ncid, jpuvar, 'long_name', 39,
     &                   'parent-to-child U-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpuvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xiu',  nf_double, 2, dimids,  xiuvar)
      ierr=nf_put_att_text(ncid, xiuvar, 'long_name', 40,
     &                 'parent-to-child U-interpolation offset X')
      ierr=nf_put_att_double(ncid, xiuvar, '_FillValue', nf_double,
     &                                                   1, -1.D0)

      ierr=nf_def_var(ncid, 'etau', nf_double, 2, dimids, etauvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     &                  'parent-to-child U-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etauvar,'_FillValue', nf_double,
     &                                                    1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'ipv', nf_int, 2, dimids, ipvvar)
      ierr=nf_put_att_text(ncid, ipvvar, 'long_name', 39,
     &                   'parent-to-child V-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jpv', nf_int, 2, dimids, jpvvar)
      ierr=nf_put_att_text(ncid, jpvvar, 'long_name', 39,
     &                   'parent-to-child V-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xiv',  nf_double, 2, dimids,  xivvar)
      ierr=nf_put_att_text(ncid, xivvar, 'long_name', 40,
     &                 'parent-to-child V-interpolation offset X')
      ierr=nf_put_att_double(ncid, xivvar, '_FillValue', nf_double,
     &                                                   1, -1.D0)

      ierr=nf_def_var(ncid, 'etav', nf_double, 2, dimids, etavvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     &                  'parent-to-child V-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavvar,'_FillValue', nf_double,
     &                                                    1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'csA', nf_double, 2, dimids, csvar)
      ierr=nf_put_att_text(ncid, csvar, 'long_name', 44,
     &          'cosine of child-parent grid angle difference')

      ierr=nf_def_var(ncid, 'snA', nf_double, 2, dimids, snvar)
      ierr=nf_put_att_text(ncid, snvar, 'long_name', 42,
     &            'sine of child-parent grid angle difference')

      ierr=nf_def_var(ncid, 'h', nf_double, 2, dimids, hvar)
      ierr=nf_put_att_text(ncid, hvar, 'long_name', 35,
     &                     'child grid native bottom topography')

      ierr=nf_def_var(ncid, 'hprnt', nf_double, 2, dimids, hpvar)
      ierr=nf_put_att_text(ncid, hpvar, 'long_name', 35,
     &                     'interpolated parent grid topography')

      ierr=nf_def_var(ncid, 'dh', nf_double, 2, dimids, dhvar)
      ierr=nf_put_att_text(ncid, dhvar, 'long_name', 34,
     &                     'child-parent topography difference')

      ierr=nf_def_var(ncid, 'kprnt', nf_real, 3, dimids, kpvar)
      ierr=nf_put_att_text(ncid, kpvar, 'long_name', 56,
     & 'index coordinate for parent-child vertical interpolation')

!>
      ierr=nf_enddef(ncid)

      ierr=nf_put_var_int(ncid, ipvar, ip)
      ierr=nf_put_var_int(ncid, jpvar, jp)
      ierr=nf_put_var_double(ncid, xivar,  xi)
      ierr=nf_put_var_double(ncid, etavar, eta)

      ierr=nf_put_var_int(ncid, ipuvar, ipu)
      ierr=nf_put_var_int(ncid, jpuvar, jpu)
      ierr=nf_put_var_double(ncid, xiuvar,  xiu)
      ierr=nf_put_var_double(ncid, etauvar, etau)

      ierr=nf_put_var_int(ncid, ipvvar, ipv)
      ierr=nf_put_var_int(ncid, jpvvar, jpv)
      ierr=nf_put_var_double(ncid, xivvar,  xiv)
      ierr=nf_put_var_double(ncid, etavvar, etav)

      ierr=nf_put_var_double(ncid, csvar, csA)
      ierr=nf_put_var_double(ncid, snvar, snA)

      ierr=nf_put_var_double(ncid, hvar,  h)
      ierr=nf_put_var_double(ncid, hpvar, hprnt)

      ierr=nf_put_var_double(ncid, kpvar, kprnt)

      allocate(dh(ncx,ncy))
      do j=1,ncy
        do i=1,ncx
          dh(i,j)=h(i,j)-hprnt(i,j)
        enddo
      enddo
      ierr=nf_put_var_double(ncid, dhvar, dh)
      deallocate(dh)
      ierr=nf_close(ncid)
      write(*,'(2x,A)') '...done'
      end
