#ifndef WRITER
# include "cppdefs.h"

      function ncdf_read(ncid, varid, record, horiz_type, A, nmax)
#else
      function ncdf_write(ncid, varid, record, horiz_type, A, nmax)
#endif

! Read/Write a floating point array from/into an input/output
! netCDF file.

! Arguments: A       real array of standard horizontal dimensions
!                                    which is to be read or written.
!            ncid    netCDF ID of in the file.
!            varid   variable ID of that variable in netCDF file.
!            record  record number.
!            type    type of the grid (RHO-, U, V, W, PSI etc.)

! Because significant portion of the code calculates starting and
! stopping indices for the subarray (which are exactly the same for
! both read and write operations) the code has been unified to
! guarantee that both functions are exactly adjoint.

      use param
      use buffer
      use grid
      use netcdf

      implicit none
#ifndef WRITER
      integer ncdf_read
#else
      integer ncdf_write
#endif
      integer ncid, varid, record, horiz_type, nmax
      real A(GLOBAL_2D_ARRAY,nmax)
      logical mask_land_data
      integer vid, i,j,k, shft, ierr
#ifdef WRITER
      integer datatype, ndims, natts, dimid(8)
      character(len=16) vname
#endif
#ifdef MASK_LAND_DATA
      real*8, parameter :: spv_set=1.D+33
#endif
#include "netcdf.inc"
#include "compute_starts_counts.h"

      if (varid > 0) then          ! Normally netCDF IDs are positive.
        vid=varid                  ! Negative "varid" is used here as
        mask_land_data=.true.      ! flag to signal that land masking
      else                         ! does not need to be applied for
        vid=-varid                 ! this variable (typically this is
        mask_land_data=.false.     ! reserved for grid variables and
      endif                        ! topography).

! Read/Write array from the disk.
!===== ===== ===== ==== === =====

#ifndef WRITER
      ierr=nf_get_vara_FTYPE(ncid, vid, start,count, buff)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,2A,3x,A,I4/)') '### ERROR: ncdf_read :: ',
     &             nf_strerror(ierr) MYID
      else
#endif

! Note that expression for "shft" is exactly the same in all five
! cases below, while application of land mask is different for the
! variables of different grid staggering; also note effectively .or.
! rather than .and. logic in setting velocity values to infinity:
! velocity components at the boundary (normal to it) are set to 0,
! while the ones fully inside (between two land points) to spv.

#ifdef MASK_LAND_DATA
        if (mask_land_data) then
          if (horiz_type == 0) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j) > 0.5) then
# ifdef WRITER
                    buff(i+shft)=A(i,j,k)
# else
                    A(i,j,k)=buff(i+shft)
# endif
                  else
# ifdef WRITER
                    buff(i+shft)=spv_set
# else
                    A(i,j,k)=0.D0
# endif
                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 1) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j)+rmask(i-1,j) > 0.5) then
# ifdef WRITER
                    buff(i+shft)=A(i,j,k)
# else
                    A(i,j,k)=buff(i+shft)
# endif
                  else
# ifdef WRITER
                    buff(i+shft)=spv_set
# else
                    A(i,j,k)=0.D0
# endif
                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 2) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j)+rmask(i,j-1) > 0.5) then
# ifdef WRITER
                    buff(i+shft)=A(i,j,k)
# else
                    A(i,j,k)=buff(i+shft)
# endif
                  else
# ifdef WRITER
                    buff(i+shft)=spv_set
# else
                    A(i,j,k)=0.D0
# endif
                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 3) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if ( rmask(i,j)+rmask(i-1,j)+rmask(i,j-1)
     &                           +rmask(i-1,j-1) > 0.5 ) then
# ifdef WRITER
                    buff(i+shft)=A(i,j,k)
# else
                    A(i,j,k)=buff(i+shft)
# endif
                  else
# ifdef WRITER
                    buff(i+shft)=spv_set
# else
                    A(i,j,k)=0.D0
# endif
                  endif
                enddo
              enddo
            enddo
          endif  !<-- horiz_type == 0,1,2,3
        else  !<-- mask_land_data
#endif
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax
# ifdef WRITER
                buff(i+shft)=A(i,j,k)
# else
                A(i,j,k)=buff(i+shft)
# endif
              enddo
            enddo
          enddo
#ifdef MASK_LAND_DATA
        endif  !<-- mask_land_data
#endif


#ifdef WRITER
        ierr=nf_put_vara_FTYPE(ncid, vid, start, count, buff)
        if (ierr /= nf90_noerr) then
          write(*,'(/1x,2A,3x,A,I4/)')  '### ERROR: ncdf_write :: ',
     &             nf_strerror(ierr) MYID
          write(*,'(12x,A,I7,3x,A,I7/12x,A,I3,7I6)') 'ncid =', ncid,
     &   'varid =', vid, 'start,count =', (start(i),count(i), i=1,4)

          i=nf_inq_var(ncid, vid,vname, datatype, ndims,dimid,natts)
          if (i == nf90_noerr) then
            write(*,'(1x,2A,1x,A,I2,2x,A,I3,2x,A,8I3)') 'vname = ',
     &             vname, 'datatype =', datatype, 'ndims =', ndims,
     &                             'dimid =', (dimid(i), i=1,ndims)
            do i=1,ndims
              j=nf_inq_dim(ncid, dimid(i), vname, k)
              if (j == nf90_noerr) write(*,'(29x,2A,I5)') vname,' =',k
            enddo
          endif
#endif
        endif

! Exchange periodic and computational margins (reader only).

#ifndef WRITER
# ifdef EXCHANGE
#  ifdef MPI
#   define EXCH_ARR_RANGE iwest,ieast,jsouth,jnorth
#  else
#   define EXCH_ARR_RANGE 1,Lm,1,Mm
#  endif
#  ifdef SOLVE3D
      call exchange_tile(EXCH_ARR_RANGE, A,nmax)
#  else
      call exchange2d_tile(EXCH_ARR_RANGE, A)
#  endif
# endif
#endif

#ifdef WRITER
      ncdf_write=ierr
#else
      ncdf_read=ierr
#endif
      end

#ifndef WRITER
# define WRITER
# include "ncdf_read_write.F"
#endif
