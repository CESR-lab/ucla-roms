! This is include file "lmd_param.h".
!----- -- ------- ---- ------------

#ifdef SOLVE3D
# ifdef LMD_MIXING
      real, parameter ::         ! Critical gradient Richardson number
     &      Ri0=0.7,                ! below which shear instabilty occurs.
#  ifdef SO_KPP
     &      nu0m=0.01,          ! Maximum viscosity and diffusivity
     &      nu0s=0.01,           ! due to shear instability [m^2/s];
#  else
     &      nu0m=50.e-4,       ! Maximum viscosity and diffusivity
     &      nu0s=50.e-4,        ! due to shear instability [m^2/s];
#  endif
     &      nuwm=1.0e-4,      ! Interior viscosity and diffusivity
     &      nuws=0.1e-4,       ! due to wave breaking, [m^2/s]

     &      lmd_nu=1.5e-6     ! Molecular viscosity [m^2/s];

#  ifdef LMD_DDMIX
      real, parameter ::         ! Value of double-diffusive density
     &      lmd_Rrho0=1.9,   ! ratio where diffusivity goes to zero
                                           ! in salt fingering.
     &      lmd_nuf=10.0e-4, ! Scaling factors for double diffusion
     &      lmd_fdd=0.7,        ! coefficient in salt fingering.
     &      lmd_tdd1=0.909,
     &      lmd_tdd2=4.6,      ! Double diffusion constants for
     &      lmd_tdd3=0.54,    ! temperature (Marmorino and Caldwell,
     &      lmd_sdd1=0.15,   ! 1976) and salinity (Fedorov, 1988)
     &      lmd_sdd2=1.85,
     &      lmd_sdd3=0.85,
     &      Smean=35.0
#  endif
# endif

# ifdef LMD_CONVEC
      real, parameter ::     
#  ifdef SO_KPP
     &   nu0c=0.01               ! convective adjustment for viscosity and
                                           ! diffusivity [m^2/s]
#  else
     &   nu0c=0.1                 ! convective adjustment for viscosity and
                                           ! diffusivity [m^2/s]
#  endif
# endif

# if defined LMD_KPP || defined LMD_BKPP
      real, parameter ::
#  ifdef SO_KPP
     &   Ricr=0.3,                 ! Critical bulk Richardson number (0.3)
#  else
     &   Ricr=0.15,               ! Critical bulk Richardson number (0.3)
                                  ! Critical bulk Richardson number. (must be decreased
                                  ! in case of a diurnal cycle, see McWilliams et al., JPO 2009)
                                  ! (recommendation 0.45 , 0.15)
#  endif
     &   Ri_inv=1./Ricr,
     &   epssfc=0.1,             ! nondimensional extent of the surface layer
     &   betaT=-0.2,             ! ratio of entrainment flux to surface buoyancy
                                           ! forcing flux (the "20% convectin rule")
     &   nubl=0.01,               ! maximum allowed boundary layer
     &   Cv=1.8,                   ! ratio of interior Brunt-Vaisala frequency
                                           ! "N" to that at the entrainment depth "he".
     &   C_MO=1.,                ! constant for computaion Monin-Obukhov depth.
#  ifdef SO_KPP
     &   C_Ek=211.,              ! constant for computating stabilization term
                                            ! due to Coriolis force (Ekman depth limit)
#  else
     &   C_Ek=258.,              ! constant for computating stabilization term
                                            ! due to Coriolis force (Ekman depth limit)
#  endif
#  if defined LMD_LIMIT_STABLE && defined SO_KPP
     &   dh0=15.,                   ! dissipation length scale [m] (Markus 1999)
     &   dc0=1.,                     ! constant (Markus 1999)
#  endif
     &   Cstar=10.,                ! proportionality coefficient parameterizing
                                            ! nonlocal transport
     &   zeta_m=-0.2,            ! Maximum stability parameters "zeta"
     &   a_m=1.257,              ! value of the 1/3 power law regime of
     &   c_m=8.360,              ! flux profile for momentum and tracers
     &   zeta_s=-1.0,             ! and coefficients of flux profile for
     &   a_s=-28.86,              ! momentum and tracers in their 1/3-power
     &   c_s=98.96                ! law regime;

      real, parameter :: r2=0.5, r3=1./3., r4=0.25
# endif
#endif
