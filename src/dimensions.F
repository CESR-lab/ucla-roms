      module dimensions

      ! ========================================
      ! Contains global and subdomain dimensions
      ! ========================================

      ! coded by Jeroen Molemaker
      use param
!      use grid

      implicit none
      private

#include "cppdefs.opt"

      integer, parameter, public :: gnx = LLm
      integer, parameter, public :: gny = MMm
      integer, parameter, public ::  nz = N

      integer, public            :: nx   ! Can vary between subdomains so not compile time parameters
      integer, public            :: ny

      ! DevinD: store global index ranges for each node for bry_extract_child.F
      integer, public :: glob_istr, glob_iend, glob_jstr, glob_jend

      public :: set_global_ij_SW_NE_corners
      public :: init_dimensions


      contains

! ----------------------------------------------------------------------
      subroutine set_global_ij_SW_NE_corners   ![ set the ij values for the south-west corner
                                               !  of subdomain based on global index numbering
      implicit none

      ! This routine includes start and end of boundary edges to include ghost point!
      ! E.g. glob_istr=0 for western edge and glob_iend=LLm+1 for eastern edge!

      if     (WESTERN_MPI_EDGE) then
        glob_istr=iSW_corn+iwest-1             ! -1 as global rho-points start from 0
      else
        glob_istr=iSW_corn+iwest
      endif

      if     (EASTERN_MPI_EDGE) then
        glob_iend=glob_istr+ieast              ! Extra end rho-point of LLm+1 is caught since no -1 unlike below
      elseif (WESTERN_MPI_EDGE) then           ! Account for fact that length is not caught by ieast fpr west
        glob_iend=glob_istr+ieast-(iwest-1)    ! boundary as can start > 1. Also, iwest-1=off_xi/2
      else
        glob_iend=glob_istr+ieast-1            ! -1 as closed interval [istr, iend], the number of
      endif                                    ! integers is ieast=iend-istr+1. Hence iend=istr+ieast-1

      if     (SOUTHERN_MPI_EDGE) then
        glob_jstr=jSW_corn+jsouth-1            ! -1 as global rho-points start from 0
      else
        glob_jstr=jSW_corn+jsouth
      endif

      if     (NORTHERN_MPI_EDGE) then
        glob_jend=glob_jstr+jnorth
      elseif (SOUTHERN_MPI_EDGE) then          ! Account for fact that length is not caught
        glob_jend=glob_jstr+jnorth-(jsouth-1)  ! by jnorth for south boundary as can start j > 1
      else                                     ! jsouth-1=off_eta/2
        glob_jend=glob_jstr+jnorth-1
      endif

!      print *, mynode,glob_istr,glob_iend,glob_jstr,glob_jend

      end subroutine set_global_ij_SW_NE_corners  !]

! ----------------------------------------------------------------------
      subroutine init_dimensions  ![
      implicit none
      ! need to do this here because Lm, Mm change in mpi_setup and are not known at compile time.

      nx = Lm
      ny = Mm

      end subroutine init_dimensions  !]

! ----------------------------------------------------------------------

      end module dimensions
