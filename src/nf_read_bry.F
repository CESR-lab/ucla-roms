#include "cppdefs.h"

      function read_bry_EW(ncid, varid, record, horiz_type, A, nmax)

! Read side boundary from an input netCDF file.
!
! Arguments:  A       real array of standard horizontal dimensions
!                             which is to be read or written.
!            ncid    netCDF ID of in the file.
!           varid    variable ID of that variable in netCDF file.
!            record  record number.
!        horiz_type   RHO-, U, V, W, PSI etc.)

      use param
      implicit none
      integer read_bry_EW, ncid, varid, record, horiz_type, nmax
      real A(0:Mm+1,nmax)
      integer i,j,k, shft, ierr, lstr, lenstr
      character(len=16) vname
#include "buffer.h"
#include "netcdf.inc"

#include "compute_starts_counts.h"

#ifdef VERBOSE
      write(*,*) 'enter read_bry_EW...' MYID
#endif

      start(1)=start(2)    ! Suppress the fist horizontal dimension
      count(1)=count(2)    ! XI: after this step the updated start(1)
      start(2)=start(3)    ! count(1) becomes corresponds to ETA-
      count(2)=count(3)    ! netCDF dimension; index 2 to vertical
      start(3)=start(4)    ! (if the variable has it) or time record;
      count(3)=count(4)    ! index 3 to time record.

! Read array from the disk.
!===== ===== ==== === =====

      ierr=nf_get_vara_FTYPE (ncid, varid, start, count, buff)
      if (ierr == nf_noerr) then
        do k=1,nmax
          shft=1-jmin + (k-1)*count(1)
          do j=jmin,jmax
            A(j,k)=buff(j+shft)
          enddo
        enddo
#ifdef MPI
        call exch_bry_EW_tile(jsouth,jnorth, A, nmax)
#elif defined NS_PERIODIC
        call exch_bry_EW_tile(     1,   Mm,  A, nmax)
#endif
      else
        ierr=nf_inq_varname(ncid, varid, vname)
        if (ierr == nf_noerr) then
          lstr=lenstr(vname)
          write(*,'(1x,2A,I6,1x,3A,2x,A,2x,A,I4)')  '### ERROR: ',
     &            'read_bry_EW :: Cannot read rec =',   record,
     &            'of variable ''',  vname(1:lstr),  '''.',
     &             nf_strerror(ierr) MYID
        else
          write(*,'(1x,2A,I6,1x,A,I4,2x,A,2x,A,I4)') '### ERROR: ',
     &            'read_bry_EW :: Cannot read rec =',   record,
     &            'of variable with netCDF ID =',  varid,
     &             nf_strerror(ierr) MYID
        endif
      endif

#ifdef VERBOSE
      write(*,*) 'leaving read_bry_EW'  MYID
#endif
      read_bry_EW=ierr
      end


      function read_bry_NS(ncid, varid, record, horiz_type, A, nmax)
      use param
      implicit none
      integer read_bry_NS, ncid, varid, record, horiz_type, nmax
      real A(0:Lm+1,nmax)
      integer  i, k, shft, ierr, lstr, lenstr
      character(len=16) vname
#include "buffer.h"
#include "netcdf.inc"

#include "compute_starts_counts.h"

#ifdef VERBOSE
      write(*,*) 'enter read_bry_NS...'   MYID
#endif
      start(2)=start(3)      ! Suppress netCDF dimension ETA. After
      count(2)=count(3)      ! this step start,count(1) correspond to
      start(3)=start(4)      ! XI-dimensions (as before); index (2) to
      count(3)=count(4)      ! vertical dimension (if present) or to
      count(3)=count(4)      ! time record; (3) (if any) to time
                             ! record.
! Read array from the disk.
!===== ===== ==== === =====

      ierr=nf_get_vara_FTYPE(ncid, varid, start, count, buff)
      if (ierr == nf_noerr) then
        do k=1,nmax
          shft=1-imin + (k-1)*count(1)
          do i=imin,imax
            A(i,k)=buff(i+shft)
          enddo
        enddo
#ifdef MPI
        call exch_bry_NS_tile(iwest,ieast, A, nmax)
#elif defined EW_PERIODIC
        call exch_bry_NS_tile(    1,   Lm, A, nmax)
#endif
      else
        ierr=nf_inq_varname (ncid, varid, vname)
        if (ierr == nf_noerr) then
          lstr=lenstr(vname)
          write(*,'(1x,2A,I6,1x,3A,2x,A,2x,A,I4)')  '### ERROR: ',
     &            'read_bry_NS :: Cannot read rec =',   record,
     &            'of variable ''',  vname(1:lstr),  '''.',
     &             nf_strerror(ierr) MYID
        else
          write(*,'(1x,2A,I6,1x,A,I4,2x,A,2x,A,I4)') '### ERROR: ',
     &            'read_bry_NS :: Cannot read rec =',   record,
     &            'of variable with netCDF ID =',  varid,
     &             nf_strerror(ierr) MYID
        endif
      endif

#ifdef VERBOSE
      write(*,*) 'leaving read_bry_NS' MYID
#endif
      read_bry_NS=ierr
      end


! The following two routines exchange ghost points of boundary arrays
! in the case of periodicity in the transversal direction or exchanging
! MPI margins. Their intended use is in the context of reading boundary
! data from a file (routines above), or creating it by other means
! (analytical). This means that the routines below can be called both
! for the whole domain (or MPI subdomain) as well as for a tile.



#if defined NS_PERIODIC || defined MPI

      subroutine exch_bry_EW_tile (jstr,jend, A, nmax)
      use param
      implicit none
      integer jstr,jend, nmax, k
      real A(0:Mm+1,nmax)
# include "buffer.h"
# if defined MPI
#  include "mpif.h"
#  include "hidden_mpi_vars.h"
      integer, parameter :: ibss=100000
      integer req(2), status(MPI_STATUS_SIZE,2), ierr
      real, dimension(N+1) :: bffr_rsv_S, bffr_snd_S,
     &                        bffr_rsv_N, bffr_snd_N
# endif

# ifdef NS_PERIODIC
#  ifdef MPI
      if (NP_ETA == 1) then
#  endif
        if (jstr == 1) then           ! Note that CPP-macros
          do k=1,nmax                 ! SOUTHERN_EDGE and NOTHERN_EDGE
            A(Mm+1,k)=A(1,k)          ! are undefined in the case of
          enddo                       ! periodicity. Also because this
        endif                         ! part is only applicable when
        if (jend == Mm) then          ! NP_ETA=1 there is no need to
          do k=1,nmax                 ! express the indices in terms
            A(0,k)=A(Mm,k)            ! of jsouth,jnorth in MPI case.
          enddo
        endif
#  ifdef MPI
      else
#  endif
# endif
# ifdef MPI
#  ifdef VERBOSE
        write(*,*) 'starting exch_bry_EW mpi exchange' MYID
#  endif
        if (SOUTH_MSG_EXCH) then
          call MPI_Irecv (bffr_rsv_S, nmax, MPI_DOUBLE_PRECISION,
     &               p_S, ibss+1, ocean_grid_comm, req(1), ierr)
        endif
        if (NORTH_MSG_EXCH) then
          call MPI_Irecv (bffr_rsv_N, nmax, MPI_DOUBLE_PRECISION,
     &               p_N, ibss+2, ocean_grid_comm, req(2), ierr)
        endif
        if (SOUTH_MSG_EXCH) then
          do k=1,nmax
            bffr_snd_S(k)=A(jstr,k)
          enddo
          call MPI_Send  (bffr_snd_S, nmax, MPI_DOUBLE_PRECISION,
     &                    p_S, ibss+2, ocean_grid_comm,    ierr)
        endif
        if (NORTH_MSG_EXCH) then
          do k=1,nmax
            bffr_snd_N(k)=A(jend,k)
          enddo
          call MPI_Send (bffr_snd_N, nmax, MPI_DOUBLE_PRECISION,
     &                   p_N, ibss+1, ocean_grid_comm,    ierr)
        endif
        if (SOUTH_MSG_EXCH) then
          call MPI_Wait (req(1),status(1,1),ierr)
          do k=1,nmax
            A(jstr-1,k)=bffr_rsv_S(k)
          enddo
        endif
        if (NORTH_MSG_EXCH) then
          call MPI_Wait (req(2),status(1,2),ierr)
          do k=1,nmax
            A(jend+1,k)=bffr_rsv_N(k)
          enddo
        endif
#  ifdef NS_PERIODIC
      endif ! <-- NP_ETA==1
#  endif
# endif
      end
#endif

#if defined EW_PERIODIC || defined MPI

      subroutine exch_bry_NS_tile (istr,iend, A, nmax)
      use param
      implicit none
      integer istr,iend, nmax, k
      real A(0:Lm+1,nmax)
# if defined MPI
#  include "mpif.h"
#  include "hidden_mpi_vars.h"
      integer, parameter :: ibss=100000
      integer req(2), status(MPI_STATUS_SIZE,2), ierr
      real, dimension(N+1) :: bffr_rsv_W, bffr_rsv_E,
     &                        bffr_snd_W, bffr_snd_E
# endif

# ifdef EW_PERIODIC
#  ifdef MPI
      if (NP_XI == 1) then
#  endif
        if (istr == 1) then           ! CPP-macros WESTERN_EDGE and
          do k=1,nmax                 ! EASTERN_EDGE are undefined in
            A(Lm+1,k)=A(1,k)          ! the case of periodicity.
          enddo                       ! Because this part is only
        endif                         ! applicable when NP_XI=1 there
        if (iend == Lm) then          ! is no need to express the
          do k=1,nmax                 ! indices in terms of iwest and
            A(0,k)=A(Lm,k)            ! ieast in MPI case.
          enddo
        endif
#  ifdef MPI
      else
#  endif
# endif
# ifdef MPI
#  ifdef VERBOSE
        write(*,*) 'starting exch_bry_NS mpi exchange' MYID
#  endif
        if (WEST_MSG_EXCH) then
          call MPI_Irecv (bffr_rsv_W, nmax, MPI_DOUBLE_PRECISION,
     &               p_W, ibss+3, ocean_grid_comm, req(1), ierr)
        endif
        if (EAST_MSG_EXCH) then
          call MPI_Irecv (bffr_rsv_E, nmax, MPI_DOUBLE_PRECISION,
     &               p_E, ibss+4, ocean_grid_comm, req(2), ierr)
        endif

        if (WEST_MSG_EXCH) then
          do k=1,nmax
            bffr_snd_W(k)=A(istr,k)
          enddo
          call MPI_Send (bffr_snd_W, nmax, MPI_DOUBLE_PRECISION,
     &                   p_W, ibss+4, ocean_grid_comm,    ierr)
        endif
        if (EAST_MSG_EXCH) then
          do k=1,nmax
            bffr_snd_E(k)=A(iend,k)
          enddo
          call MPI_Send (bffr_snd_E, nmax, MPI_DOUBLE_PRECISION,
     &                   p_E, ibss+3, ocean_grid_comm,    ierr)
        endif

        if (WEST_MSG_EXCH) then
          call MPI_Wait (req(1),status(1,1),ierr)
          do k=1,nmax
            A(istr-1,k)=bffr_rsv_W(k)
          enddo
        endif
        if (EAST_MSG_EXCH) then
          call MPI_Wait (req(2),status(1,2),ierr)
          do k=1,nmax
            A(iend+1,k)=bffr_rsv_E(k)
          enddo
        endif
#  ifdef EW_PERIODIC
      endif !<-- NP_XI==1
#  endif
# endif
      end
#endif
