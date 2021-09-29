/* Fill two rows of ghost points associated with periodic boundaries,
 if any, and if there is no MPI partitioning in the periodic direction,
 then call MPI halo exchange routines to fill up computational margins.
 This can be done for one, two, three, or four arrays simultaneously,
 so the actual messages in MPI exchange routines are formed packing
 data from different arrays resulting in larger, but fewer messages in
 comparison with communicating each array separately. For convenience
 there are two sets of routines named

        exchange2d_tile   exchange2d_X_tile
 and
        exchange_tile     exchange_X_tile

 where X=2,3,4 (number of arrays, if more than one) with "2d" is
 strictly for 2D arrays, while the second set can communicate both 2D
 and 3D arrays in any combination by explicitly specifying depth
 dimension after each array argument. Note that neither these routines,
 nor the MPI exchange calls are designed to distinguish between RHO-,
 U-, V- or PSI-type variables. CPP-switch BAKWD_COMPAT is provided for
 backward compatibility with the previous versions of the code, where
 the array types were distingiushed and the arrays have to be exchanged
 one at at time. Normally this switch should be undefined.   */

#ifndef NARGS
# include "cppdefs.opt"
!--># define BAKWD_COMPAT


/*
# define mpi_exchange_tile   mpi_exchange4_tile
# define mpi_exchange_2_tile mpi_exchange4_2_tile
# define mpi_exchange_3_tile mpi_exchange4_3_tile
# define mpi_exchange_4_tile mpi_exchange4_4_tile
*/

# define mpi_exchange_tile   mpi_exchange8_tile
# define mpi_exchange_2_tile mpi_exchange8_2_tile
# define mpi_exchange_3_tile mpi_exchange8_3_tile
# define mpi_exchange_4_tile mpi_exchange8_4_tile



#endif
#ifdef EXCHANGE
# ifndef NARGS
      subroutine exchange2d_tile (istr,iend,jstr,jend,  A)
# elif NARGS == 2
      subroutine exchange2d_2_tile (istr,iend,jstr,jend, A,B)
# elif NARGS == 3
      subroutine exchange2d_3_tile (istr,iend,jstr,jend, A,B,C)
# elif NARGS == 4
      subroutine exchange2d_4_tile (istr,iend,jstr,jend, A,B,C,D)
# endif
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY)
# if NARGS > 1
      real B(GLOBAL_2D_ARRAY)
# endif
# if NARGS > 2
      real C(GLOBAL_2D_ARRAY)
# endif
# if NARGS > 3
      real D(GLOBAL_2D_ARRAY)
# endif
# ifdef EW_PERIODIC
      integer j
# endif
# ifdef NS_PERIODIC
      integer i
# endif

# include "compute_auxiliary_bounds.h"

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE jstr,jend
#  else
#   define J_RANGE jstrR,jendR
#  endif
#  ifdef MPI
      if (NP_XI == 1) then            !!!!!!  WARNING !!!!!!!
#  endif
        if (istr == 1) then           ! It should be noted that the
          do j=J_RANGE                ! four CPP-macros,
            A(Lm+1,j)=A(1,j)          ! WESTERN_EDGE,  EASTERN_EDGE,
            A(Lm+2,j)=A(2,j)          ! SOUTHERN_EDGE, NORTHERN_EDGE,
#  if NARGS > 1
            B(Lm+1,j)=B(1,j)          ! in "set_global_definitions.h"
            B(Lm+2,j)=B(2,j)          ! are specifically designed to
#  endif
#  if NARGS > 2
            C(Lm+1,j)=C(1,j)          ! be undefined in the case of
            C(Lm+2,j)=C(2,j)          ! periodicity in the respective
#  endif
#  if NARGS > 3
            D(Lm+1,j)=D(1,j)          ! horizontal dimension, so the
            D(Lm+2,j)=D(2,j)          ! use of these macros is avoided
#  endif
          enddo                       ! throughout this file.
        endif

        if (iend == Lm) then
          do j=J_RANGE
            A(-1,j)=A(Lm-1,j)
            A( 0,j)=A(Lm  ,j)
#  if NARGS > 1
            B(-1,j)=B(Lm-1,j)
            B( 0,j)=B(Lm  ,j)
#  endif
#  if NARGS > 2
            C(-1,j)=C(Lm-1,j)
            C( 0,j)=C(Lm  ,j)
#  endif
#  if NARGS > 3
            D(-1,j)=D(Lm-1,j)
            D( 0,j)=D(Lm  ,j)
#  endif
          enddo
        endif
#  ifdef MPI
      endif
#  endif
#  undef J_RANGE
# endif            /* EW_PERIODIC */
 
# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE istr,iend
#  else
#   define I_RANGE istrR,iendR
#  endif
#  ifdef MPI
      if (NP_ETA == 1) then
#  endif
        if (jstr == 1) then
          do i=I_RANGE
            A(i,Mm+1)=A(i,1)
            A(i,Mm+2)=A(i,2)
#  if NARGS > 1
            B(i,Mm+1)=B(i,1)
            B(i,Mm+2)=B(i,2)
#  endif
#  if NARGS > 2
            C(i,Mm+1)=C(i,1)
            C(i,Mm+2)=C(i,2)
#  endif
#  if NARGS > 3
            D(i,Mm+1)=D(i,1)
            D(i,Mm+2)=D(i,2)
#  endif
          enddo
        endif

        if (jend == Mm) then
          do i=I_RANGE
            A(i,-1)=A(i,Mm-1)
            A(i, 0)=A(i,Mm  )
#  if NARGS > 1
            B(i,-1)=B(i,Mm-1)
            B(i, 0)=B(i,Mm  )
#  endif
#  if NARGS > 2
            C(i,-1)=C(i,Mm-1)
            C(i, 0)=C(i,Mm  )
#  endif
#  if NARGS > 3
            D(i,-1)=D(i,Mm-1)
            D(i, 0)=D(i,Mm  )
#  endif
          enddo
        endif
#  ifdef MPI
      endif
#  endif
#  undef I_RANGE
# endif               /* NS_PERIODIC */
 
# if defined EW_PERIODIC && defined NS_PERIODIC
#  ifdef MPI
      if (NP_XI == 1 .and. NP_ETA == 1) then
#  endif
        if (istr == 1 .and. jstr == 1) then
          A(Lm+1,Mm+1)=A(1,1) ; A(Lm+1,Mm+2)=A(1,2)
          A(Lm+2,Mm+1)=A(2,1) ; A(Lm+2,Mm+2)=A(2,2)
#  if NARGS > 1
          B(Lm+1,Mm+1)=B(1,1) ; B(Lm+1,Mm+2)=B(1,2)
          B(Lm+2,Mm+1)=B(2,1) ; B(Lm+2,Mm+2)=B(2,2)
#  endif
#  if NARGS > 2
          C(Lm+1,Mm+1)=C(1,1) ; C(Lm+1,Mm+2)=C(1,2)
          C(Lm+2,Mm+1)=C(2,1) ; C(Lm+2,Mm+2)=C(2,2)
#  endif
#  if NARGS > 3
          D(Lm+1,Mm+1)=D(1,1) ; D(Lm+1,Mm+2)=D(1,2)
          D(Lm+2,Mm+1)=D(2,1) ; D(Lm+2,Mm+2)=D(2,2)
#  endif
        endif

        if (iend == Lm .and. jstr == 1) then
          A(-1,Mm+1)=A(Lm-1,1) ; A(0,Mm+1)=A(Lm,1)
          A(-1,Mm+2)=A(Lm-1,2) ; A(0,Mm+2)=A(Lm,2)
#  if NARGS > 1
          B(-1,Mm+1)=B(Lm-1,1) ; B(0,Mm+1)=B(Lm,1)
          B(-1,Mm+2)=B(Lm-1,2) ; B(0,Mm+2)=B(Lm,2)
#  endif
#  if NARGS > 2
          C(-1,Mm+1)=C(Lm-1,1) ; C(0,Mm+1)=C(Lm,1)
          C(-1,Mm+2)=C(Lm-1,2) ; C(0,Mm+2)=C(Lm,2)
#  endif
#  if NARGS > 3
          D(-1,Mm+1)=D(Lm-1,1) ; D(0,Mm+1)=D(Lm,1)
          D(-1,Mm+2)=D(Lm-1,2) ; D(0,Mm+2)=D(Lm,2)
#  endif
        endif

        if (istr == 1 .and. jend == Mm) then
          A(Lm+1,-1)=A(1,Mm-1) ; A(Lm+1,0)=A(1,Mm)
          A(Lm+2,-1)=A(2,Mm-1) ; A(Lm+2,0)=A(2,Mm)
#  if NARGS > 1
          B(Lm+1,-1)=B(1,Mm-1) ; B(Lm+1,0)=B(1,Mm)
          B(Lm+2,-1)=B(2,Mm-1) ; B(Lm+2,0)=B(2,Mm)
#  endif
#  if NARGS > 2
          C(Lm+1,-1)=C(1,Mm-1) ; C(Lm+1,0)=C(1,Mm)
          C(Lm+2,-1)=C(2,Mm-1) ; C(Lm+2,0)=C(2,Mm)
#  endif
#  if NARGS > 3
          D(Lm+1,-1)=D(1,Mm-1) ; D(Lm+1,0)=D(1,Mm)
          D(Lm+2,-1)=D(2,Mm-1) ; D(Lm+2,0)=D(2,Mm)
#  endif
        endif

        if (iend == Lm .and. jend == Mm) then
          A(-1,-1)=A(Lm-1,Mm-1) ; A(0,-1)=A(Lm,Mm-1)
          A(-1, 0)=A(Lm-1,Mm  ) ; A(0, 0)=A(Lm,Mm  )
#  if NARGS > 1
          B(-1,-1)=B(Lm-1,Mm-1) ; B(0,-1)=B(Lm,Mm-1)
          B(-1, 0)=B(Lm-1,Mm  ) ; B(0, 0)=B(Lm,Mm  )
#  endif
#  if NARGS > 2
          C(-1,-1)=C(Lm-1,Mm-1) ; C(0,-1)=C(Lm,Mm-1)
          C(-1, 0)=C(Lm-1,Mm  ) ; C(0, 0)=C(Lm,Mm  )
#  endif
#  if NARGS > 3
          D(-1,-1)=D(Lm-1,Mm-1) ; D(0,-1)=D(Lm,Mm-1)
          D(-1, 0)=D(Lm-1,Mm  ) ; D(0, 0)=D(Lm,Mm  )
#  endif
        endif
#  ifdef MPI
      endif
#  endif
# endif
# ifdef MPI
#  ifndef NARGS
      call mpi_exchange_tile(istr,iend,jstr,jend, A,1)
#  elif NARGS == 2
      call mpi_exchange_2_tile(istr,iend,jstr,jend, A,1, B,1)
#  elif NARGS == 3
      call mpi_exchange_3_tile(istr,iend,jstr,jend, A,1, B,1, C,1)
#  elif NARGS == 4
      call mpi_exchange_4_tile(istr,iend,jstr,jend, A,1, B,1,
     &                                              C,1, D,1)
#  endif
# endif
      end


# ifdef SOLVE3D
#  ifndef NARGS
      subroutine exchange_tile(istr,iend,jstr,jend,  A, nmaxA)
#  elif NARGS == 2
      subroutine exchange_2_tile(istr,iend,jstr,jend, A, nmaxA,
     &                                                B, nmaxB)
#  elif NARGS == 3
      subroutine exchange_3_tile(istr,iend,jstr,jend, A, nmaxA,
     &                                      B, nmaxB, C, nmaxC)
#  elif NARGS == 4
      subroutine exchange_4_tile(istr,iend,jstr,jend, A, nmaxA,
     &                            B, nmaxB, C, nmaxC, D, nmaxD)
#  endif
      use param
      implicit none
      integer istr,iend,jstr,jend, nmaxA
      real A(GLOBAL_2D_ARRAY,nmaxA)
#  if NARGS > 1
      integer nmaxB
      real B(GLOBAL_2D_ARRAY,nmaxB)
#  endif
#  if NARGS > 2
      integer nmaxC
      real C(GLOBAL_2D_ARRAY,nmaxC)
#  endif
#  if NARGS > 3
      integer nmaxD
      real D(GLOBAL_2D_ARRAY,nmaxD)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
      integer k
#   ifdef EW_PERIODIC
     &       , j 
#   endif
#   ifdef NS_PERIODIC
     &       , i
#   endif
#  endif

#  include "compute_auxiliary_bounds.h"

#  ifdef EW_PERIODIC
#   ifdef NS_PERIODIC
#    define J_RANGE jstr,jend
#   else
#    define J_RANGE jstrR,jendR
#   endif
#   ifdef MPI
      if (NP_XI == 1) then
#   endif
        if (istr == 1) then
          do k=1,nmaxA
            do j=J_RANGE
              A(Lm+1,j,k)=A(1,j,k)
              A(Lm+2,j,k)=A(2,j,k)
            enddo
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            do j=J_RANGE
              B(Lm+1,j,k)=B(1,j,k)
              B(Lm+2,j,k)=B(2,j,k)
            enddo
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            do j=J_RANGE
              C(Lm+1,j,k)=C(1,j,k)
              C(Lm+2,j,k)=C(2,j,k)
            enddo
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            do j=J_RANGE
              D(Lm+1,j,k)=D(1,j,k)
              D(Lm+2,j,k)=D(2,j,k)
            enddo
          enddo
#   endif
        endif

        if (iend == Lm) then
          do k=1,nmaxA
            do j=J_RANGE
              A(-1,j,k)=A(Lm-1,j,k)
              A( 0,j,k)=A(Lm  ,j,k)
            enddo
          enddo

#   if NARGS > 1
          do k=1,nmaxB
            do j=J_RANGE
              B(-1,j,k)=B(Lm-1,j,k)
              B( 0,j,k)=B(Lm  ,j,k)
            enddo
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            do j=J_RANGE
              C(-1,j,k)=C(Lm-1,j,k)
              C( 0,j,k)=C(Lm  ,j,k)
            enddo
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            do j=J_RANGE
              D(-1,j,k)=D(Lm-1,j,k)
              D( 0,j,k)=D(Lm  ,j,k)
            enddo
          enddo
#   endif
        endif
#   ifdef MPI
      endif
#   endif
#   undef J_RANGE
#  endif            /* EW_PERIODIC */
 
#  ifdef NS_PERIODIC
#   ifdef EW_PERIODIC
#    define I_RANGE istr,iend
#   else
#    define I_RANGE istrR,iendR
#   endif
#   ifdef MPI
      if (NP_ETA == 1) then
#   endif
        if (jstr == 1) then
          do k=1,nmaxA
            do i=I_RANGE
              A(i,Mm+1,k)=A(i,1,k)
              A(i,Mm+2,k)=A(i,2,k)
            enddo
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            do i=I_RANGE
              B(i,Mm+1,k)=B(i,1,k)
              B(i,Mm+2,k)=B(i,2,k)
            enddo
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            do i=I_RANGE
              C(i,Mm+1,k)=C(i,1,k)
              C(i,Mm+2,k)=C(i,2,k)
            enddo
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            do i=I_RANGE
              D(i,Mm+1,k)=D(i,1,k)
              D(i,Mm+2,k)=D(i,2,k)
            enddo
          enddo
#   endif
        endif

        if (jend == Mm) then
          do k=1,nmaxA
            do i=I_RANGE
              A(i,-1,k)=A(i,Mm-1,k)
              A(i, 0,k)=A(i,Mm  ,k)
            enddo
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            do i=I_RANGE
              B(i,-1,k)=B(i,Mm-1,k)
              B(i, 0,k)=B(i,Mm  ,k)
            enddo
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            do i=I_RANGE
              C(i,-1,k)=C(i,Mm-1,k)
              C(i, 0,k)=C(i,Mm  ,k)
            enddo
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            do i=I_RANGE
              D(i,-1,k)=D(i,Mm-1,k)
              D(i, 0,k)=D(i,Mm  ,k)
            enddo
          enddo
#   endif
        endif
#   ifdef MPI
      endif
#   endif
#   undef I_RANGE
#  endif               /* NS_PERIODIC */
 
#  if defined EW_PERIODIC && defined NS_PERIODIC
#   ifdef MPI
      if (NP_XI == 1 .and. NP_ETA == 1) then
#   endif
        if (istr == 1 .and. jstr == 1) then
          do k=1,nmaxA
            A(Lm+1,Mm+1,k)=A(1,1,k) ; A(Lm+1,Mm+2,k)=A(1,2,k)
            A(Lm+2,Mm+1,k)=A(2,1,k) ; A(Lm+2,Mm+2,k)=A(2,2,k)
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            B(Lm+1,Mm+1,k)=B(1,1,k) ; B(Lm+1,Mm+2,k)=B(1,2,k)
            B(Lm+2,Mm+1,k)=B(2,1,k) ; B(Lm+2,Mm+2,k)=B(2,2,k)
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            C(Lm+1,Mm+1,k)=C(1,1,k) ; C(Lm+1,Mm+2,k)=C(1,2,k)
            C(Lm+2,Mm+1,k)=C(2,1,k) ; C(Lm+2,Mm+2,k)=C(2,2,k)
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            D(Lm+1,Mm+1,k)=D(1,1,k) ; D(Lm+1,Mm+2,k)=D(1,2,k)
            D(Lm+2,Mm+1,k)=D(2,1,k) ; D(Lm+2,Mm+2,k)=D(2,2,k)
          enddo
#   endif
        endif

        if (iend == Lm .and. jstr == 1) then
          do k=1,nmaxA
            A(-1,Mm+1,k)=A(Lm-1,1,k) ; A(0,Mm+1,k)=A(Lm,1,k)
            A(-1,Mm+2,k)=A(Lm-1,2,k) ; A(0,Mm+2,k)=A(Lm,2,k)
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            B(-1,Mm+1,k)=B(Lm-1,1,k) ; B(0,Mm+1,k)=B(Lm,1,k)
            B(-1,Mm+2,k)=B(Lm-1,2,k) ; B(0,Mm+2,k)=B(Lm,2,k)
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            C(-1,Mm+1,k)=C(Lm-1,1,k) ; C(0,Mm+1,k)=C(Lm,1,k)
            C(-1,Mm+2,k)=C(Lm-1,2,k) ; C(0,Mm+2,k)=C(Lm,2,k)
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            D(-1,Mm+1,k)=D(Lm-1,1,k) ; D(0,Mm+1,k)=D(Lm,1,k)
            D(-1,Mm+2,k)=D(Lm-1,2,k) ; D(0,Mm+2,k)=D(Lm,2,k)
          enddo
#   endif
        endif

        if (istr == 1 .and. jend == Mm) then
          do k=1,nmaxA
            A(Lm+1,-1,k)=A(1,Mm-1,k) ; A(Lm+1,0,k)=A(1,Mm,k)
            A(Lm+2,-1,k)=A(2,Mm-1,k) ; A(Lm+2,0,k)=A(2,Mm,k)
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            B(Lm+1,-1,k)=B(1,Mm-1,k) ; B(Lm+1,0,k)=B(1,Mm,k)
            B(Lm+2,-1,k)=B(2,Mm-1,k) ; B(Lm+2,0,k)=B(2,Mm,k)
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            C(Lm+1,-1,k)=C(1,Mm-1,k) ; C(Lm+1,0,k)=C(1,Mm,k)
            C(Lm+2,-1,k)=C(2,Mm-1,k) ; C(Lm+2,0,k)=C(2,Mm,k)
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            D(Lm+1,-1,k)=D(1,Mm-1,k) ; D(Lm+1,0,k)=D(1,Mm,k)
            D(Lm+2,-1,k)=D(2,Mm-1,k) ; D(Lm+2,0,k)=D(2,Mm,k)
          enddo
#   endif
        endif

        if (iend == Lm .and. jend == Mm) then
          do k=1,nmaxA
            A(-1,-1,k)=A(Lm-1,Mm-1,k) ; A(0,-1,k)=A(Lm,Mm-1,k)
            A(-1, 0,k)=A(Lm-1,Mm  ,k) ; A(0, 0,k)=A(Lm,Mm  ,k)
          enddo
#   if NARGS > 1
          do k=1,nmaxB
            B(-1,-1,k)=B(Lm-1,Mm-1,k) ; B(0,-1,k)=B(Lm,Mm-1,k)
            B(-1, 0,k)=B(Lm-1,Mm  ,k) ; B(0, 0,k)=B(Lm,Mm  ,k)
          enddo
#   endif
#   if NARGS > 2
          do k=1,nmaxC
            C(-1,-1,k)=C(Lm-1,Mm-1,k) ; C(0,-1,k)=C(Lm,Mm-1,k)
            C(-1, 0,k)=C(Lm-1,Mm  ,k) ; C(0, 0,k)=C(Lm,Mm  ,k)
          enddo
#   endif
#   if NARGS > 3
          do k=1,nmaxD
            D(-1,-1,k)=D(Lm-1,Mm-1,k) ; D(0,-1,k)=D(Lm,Mm-1,k)
            D(-1, 0,k)=D(Lm-1,Mm  ,k) ; D(0, 0,k)=D(Lm,Mm  ,k)
          enddo
#   endif
        endif
#   ifdef MPI
      endif
#   endif
#  endif

#  ifdef MPI
#   ifndef NARGS
      call mpi_exchange_tile (istr,iend,jstr,jend, A, nmaxA)
#   elif NARGS == 2
      call mpi_exchange_2_tile(istr,iend,jstr,jend, A, nmaxA,
     &                                              B, nmaxB)
#   elif NARGS == 3
      call mpi_exchange_3_tile(istr,iend,jstr,jend, A, nmaxA,
     &                                    B, nmaxB, C, nmaxC)
#   elif NARGS == 4
      call mpi_exchange_4_tile(istr,iend,jstr,jend, A, nmaxA,
     &                          B, nmaxB, C, nmaxC, D, nmaxD)
#   endif

/* Legacy version, one mpi_exchange call for each array:

      call mpi_exchange_tile(istr,iend,jstr,jend, A, nmaxA)
#   if NARGS > 1
      call mpi_exchange_tile(istr,iend,jstr,jend, B, nmaxB)
#   endif
#   if NARGS > 2
      call mpi_exchange_tile(istr,iend,jstr,jend, C, nmaxC)
#   endif
#   if NARGS > 3
      call mpi_exchange_tile(istr,iend,jstr,jend, D, nmaxD)
#   endif
*/
#  endif
      end
# endif /* SOLVE3D */

# ifndef NARGS
#  define NARGS 2
# elif NARGS == 2
#  undef NARGS
#  define NARGS 3
# elif NARGS == 3
#  undef NARGS
#  define NARGS 4
# elif NARGS == 4
#  undef NARGS
# endif
# ifdef NARGS
#  include "exchange.F"
# else
#  ifdef BAKWD_COMPAT

! The following nine subroutines are purely for the backward
! compatibility with the previous versions of the code merely
! to facilitate transition. Normally CPP-switch BAKWD_COMPAT
! should not be activated.

      subroutine exchange_r2d_tile(istr,iend,jstr,jend, A) 
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY)
      call exchange2d_tile(istr,iend,jstr,jend, A)
      end

      subroutine exchange_u2d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY)
      call exchange2d_tile(istr,iend,jstr,jend, A)
      end

      subroutine exchange_v2d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY)
      call exchange2d_tile(istr,iend,jstr,jend, A)
      end

      subroutine exchange_p2d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY)
      call exchange2d_tile(istr,iend,jstr,jend, A)
      end

#   ifdef SOLVE3D
      subroutine exchange_r3d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY,N)
      call exchange_tile(istr,iend,jstr,jend, A,N)
      end

      subroutine exchange_u3d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY,N)
      call exchange_tile(istr,iend,jstr,jend, A,N)
      end

      subroutine exchange_v3d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY,N)
      call exchange_tile(istr,iend,jstr,jend, A,N)
      end

      subroutine exchange_p3d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY,N)
      call exchange_tile(istr,iend,jstr,jend, A,N)
      end

      subroutine exchange_w3d_tile(istr,iend,jstr,jend, A)
      use param
      implicit none
      integer istr,iend,jstr,jend
      real A(GLOBAL_2D_ARRAY,0:N)
      call exchange_tile(istr,iend,jstr,jend, A,N+1)
      end

#   endif  /* SOLVE3D */
#  endif  /* BAKWD_COMPAT */ 
# endif  /* NARGS */
#else
      subroutine exchange_empty
      end
#endif /* EXCHANGE */
