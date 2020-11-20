! Tcoef, T0   Coefficients for linear Equation of State 
! Scoef, S0     rho = Tcoef*(T-T0) + Scoef*(S-S0)
!
#ifdef SOLVE3D
# ifndef NONLIN_EOS
      real Tcoef, T0
      common /eos_pars/ Tcoef, T0
#  ifdef SALINITY
      real Scoef, S0 
      common /eos_pars/ Scoef, S0
#  endif
# endif


# ifdef SPLIT_EOS
      real rho1(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE rho1(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /eos_rho1/ rho1

      real qp1(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE qp1(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /eos_qp1/ qp1
      real, parameter :: qp2=0.0000172
# else
      real rho(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE rho(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /eos_rho/ rho
# endif
# ifdef ADV_ISONEUTRAL
      real dRdx(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE dRdx(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real dRde(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE dRde(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_drdx/dRdx /ocean_drde/dRde
      real idRz(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE idRz(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_drdz/idRz
# endif
#endif
