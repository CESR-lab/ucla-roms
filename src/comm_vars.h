! The following common block contains a set of globally accessible
! variables needed for data exchange between OpenMP threads working
! on different subdomains.
!
! Global summation variables are declared with 16-byte precision to
! avoid accumulation of roundoff errors (subject to whether compiler
! suppors 16-byte real), since roundoff error depends on the order of
! summation, which is indeterministic in the case of summation between
! the parallel threads; not doing so would make it impossible to pass
! an ETALON CHECK test if there is a feedback of global sums into the
! dynamics of the model, such as in the case when global mass
! conservation is enforced.
!
! WARNING: FRAGILE ALIGNMENT SEQUENCE. In the following common block
! real objects are grouped in pairs and integer*4 are grouped in quads
! to guarantee that 16-Byte objects are aligned in 16-Byte boundaries
! and 8-Byte objects are aligned in 8-Byte boundaries.   Removal or
! introduction of variables with violation of these parity rules, as
! well as changing the sequence of variables in the common block may
! cause violation of alignment.

      real*QUAD area, volume, bc_crss
      common /comm_vars/ area, volume, bc_crss

      real hmin,hmax, grdmin,grdmax, rx0,rx1, Cg_min,Cg_max, Cu_Cor
      common /comm_vars/ hmin,hmax, grdmin,grdmax, rx0,rx1,
     &                                        Cg_min,Cg_max, Cu_Cor
      real*4 cpu_all(2)
      integer trd_count, tile_count, bc_count, mcheck, first_time
      common /comm_vars/ cpu_all, trd_count,
     &                   tile_count, bc_count, mcheck, first_time
