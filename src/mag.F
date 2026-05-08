      module mag
      ! -----------------------------------------------------------
      !                 ROMS module for:
      !                    MACROALGAE (=seaweed=kelp)
      !         
      !                 Originally coded following MAG model in Friedeer
      !                 et al. 2022 (FMARS), however, recently cleared
      !                 and made simpler to simulate seaweed CDR with
      !                 prescribed biomass and surface injections of
      !                 seaweed POC
      !                
      !                 Future revisions of this module will build back
      !                 in complexity, i.e., growth dynamics
      !            
      !             Coded by D. Dauhajre with help from P.
      !             Damien, J. Molemaker, and D. Dollery
      ! -----------------------------------------------------------

#include "cppdefs.opt"
! CPP flag MAG covering entire module so that the module is empty if
! MAG is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined MAG /* MAG for whole module! */


      use dimensions
      use nc_read_write
      use roms_read_write
      use netcdf
      use param
      use grid
      use ocean_vars
      use scalars
      use mixing
      use tracers
      use roms_mpi
      !use mpi_exchanges      

      implicit none


!  standard user-inputs (switches for forcing and writing)
#include "mag.opt"

!    constants that are not generally changed (e.g., drag coefficient)
#include "mag_params.opt"
      

      ! Make all variable private to this module unless public specified
      private
     
      ! save -> Make all public variables global
      ! (might not be needed depending on compiler)
      save

      ! -------------------------------------------------------------------
      !                    PRIMARY  MAG VARIABLES
      !        ***********************************************
      ! sflx_pockelp_mass  | surface flux of kelp biomass (for sinking) | [g C / m^2 / day] 
      ! algae_b2d          | vert. integrated biomass | [g-dry / m^2]
      ! algae_b            | 3D biomass               | [g-dry / m^3]
      ! algae_lad          | leaf-area density        | [1/m]
      ! algae_shape        | vert. biomass shape func.| [1/m0]
      ! algae_length       | length(height) of algae  | [m]
      ! algae_mask         | mask [1=algae, 0=none]   | no units
      ! algae_Nf           | fixed nitrogen of seaweed | [mg N / m^3]
      ! algae_npp          | vert. int npp             | [g-dry / m^2 / sec]
      ! algae_loss          | vert int losses          | [g-dry/ m^2 /sec]
      !       ************************************************
      ! ---------------------------------------------------------------------

      ! This is used to write PAR_lay into 3D variable in BEC
      ! which MAG needs to do growth/loss based on 3D PAR
      !         CAN WRITE THIS AS WELL AS A DIAGNOSTIC
            ! ALGAE AGE
      real,allocatable,dimension(:,:)          :: algae_age ![seconds]
      real,allocatable,dimension(:,:)          :: sflx_pockelp_mass
      real(8),allocatable,dimension(:,:)          :: algae_b2d
      real,allocatable,dimension(:,:,:)        :: algae_b
      real,allocatable,dimension(:,:,:),public :: algae_lad
      real,allocatable,dimension(:,:,:)        :: algae_shape
      real,allocatable,dimension(:,:)          :: algae_length
      real,allocatable,dimension(:,:)          :: algae_length_sigcorr !respects grid
      real,allocatable,dimension(:,:)          :: algae_mask
      real,allocatable,dimension(:,:,:),public :: algae_Nf !modifies PAR in bgc_ecosys_bec2.F 
      real(8),allocatable,dimension(:,:)         :: algae_npp
      real(8),allocatable,dimension(:,:)         :: algae_loss
      real(8),allocatable,dimension(:,:)         :: algae_harv




      ! MAG needs 3D PAR to diagnose growth/loss
      !BEC does not have this (it writes 3D PAR with diagnostics)
      ! so in BEC, we grab 3D PAR from PAR_lay (a 1D work array
      ! basically) 
      real,allocatable,dimension(:,:,:),public:: PAR_3d_mag


      !-----------------------------------------------------------
      !                 BEC COUPLING TERMS (bgc source/sinks from MAG)
      !            computed locally here during step_mag 
      !            *remin of sinking biomass computed in BEC (not here)
      !             as the sinking biomass is a BEC tracer (POCkelp) 
      !        ***********************************************
      ! algae_uptno3   | nitrate uptake by algae          | [mmol N / m^3 / s]     
      ! algae_exuL     | loss due to exudation of DOM     | [mmol N / m^3 / s]     
      ! algae_spomL    | loss of small particulates       | [mmol N / m^3 / s]     
      ! algae_frondL   | loss of fronds (big particulates)| [mmol N / m^3 / s]    
      ! algae_harvL    | loss of harvest (processed kelp) | [mmol N / m^3 / s]    
      !       ************************************************
      ! ---------------------------------------------------------------------
      real,allocatable,dimension(:,:,:),public :: algae_uptno3 
      real,allocatable,dimension(:,:,:),public :: algae_exuL   
      real,allocatable,dimension(:,:,:),public :: algae_spomL 
      real,allocatable,dimension(:,:,:),public :: algae_frondL
      real,allocatable,dimension(:,:,:),public :: algae_harvL
      real,allocatable,dimension(:,:,:),public :: remin_lrgpocalgae
      real,allocatable,dimension(:,:,:),public :: remin_smpocalgae
      real,allocatable,dimension(:,:,:),public :: remin_hrvpocalgae
      real,allocatable,dimension(:,:,:),public :: remin_docalgae
      real,allocatable,dimension(:,:,:),public :: remin_docrefalgae




      !Averages for writing (same units as above)
      real,allocatable,dimension(:,:)          :: algae_b2d_avg
      real,allocatable,dimension(:,:,:)        :: algae_b_avg
      real,allocatable,dimension(:,:,:)        :: algae_lad_avg
      real,allocatable,dimension(:,:)           :: algae_length_avg
      real,allocatable,dimension(:,:,:)        :: algae_shape_avg
      real,allocatable,dimension(:,:,:)        :: algae_uptno3_avg
      real,allocatable,dimension(:,:,:)        :: algae_exuL_avg
      real,allocatable,dimension(:,:,:)        :: algae_spomL_avg
      real,allocatable,dimension(:,:,:)        :: algae_frondL_avg
      real,allocatable,dimension(:,:,:)        :: algae_harvL_avg
      real,allocatable,dimension(:,:,:)        :: PAR_3d_mag_avg
      real,allocatable,dimension(:,:)          :: algae_npp_avg
      real,allocatable,dimension(:,:)          :: algae_loss_avg
      real, allocatable,dimension(:,:)         :: algae_harv_avg
      real,allocatable,dimension(:,:,:)        :: remin_lrgpocalgae_avg
      real,allocatable,dimension(:,:,:)        :: remin_smpocalgae_avg
      real,allocatable,dimension(:,:,:)        :: remin_hrvpocalgae_avg
      real,allocatable,dimension(:,:,:)        :: remin_docalgae_avg
      real,allocatable,dimension(:,:,:)        :: remin_docrefalgae_avg



      !surface flux of kelp biomass in units of [mmol C / m^2 / sec]
      ! this will go into the kelpPOChrv pool of POC ('harvested') 
      real,allocatable,dimension(:,:),public :: pockelp_prod_surf !used in bgc_ecosys_bec2.F

      
      !-----------------------------------------------------------
      ! 2D ARRAYS DEFINING BOTTOM OF ALGAE (cultivation depth or
      ! sea-floor if natural)
      !------------------------------------------------------------
      real, allocatable, dimension(:,:)   :: algae_zbot  !bottom dept in [m] (POSITIVE)
      integer, allocatable,dimension(:,:) :: algae_kbot  ! vert. rho-level at bottom of plant



      !-------------------------------------
      ! Boolean for initializing kelp
      ! gets set to true in set_mag_frc after kelp initialized
      ! public b/c it gets used in step3d_t_ISO so that MAG evolves only
      ! after initialization
      !--------------------------------------
      logical, public :: kelp_initialized


      !Variables that are used locally in compute_horiz_rhs_uv_terms.h
      !for momentum sink
#ifdef ALGAE_DRAG
      ! Volume normalized drag term that is subtracted from ru (u-pts)
      real, allocatable, dimension (:,:,:),public :: vegdrag_u
      !Volume normalized drag term that is subtracted from rv (v-pts)
      real, allocatable, dimension (:,:,:),public :: vegdrag_v
      !Scratch variables to make code easier to read
      !         Eventually get rid of these for efficiency
      real, public :: vegcoeffu
      real, public  :: vegcoeffv
      real, public  :: vegvolu
      real, public  :: vegvolv
      real, public  :: v_ri1
      real, public  :: v_ri2
      real, public  :: u_rj1
      real, public :: u_rj2
      real, public :: velocmagu
      real, public :: velocmagv
#endif



      ! -------------------------------------------------------------------
      !
      !                      DECALRE THINGS FOR NETCDF OUTPUT
      !
      ! ---------------------------------------------------------------------
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Variable names and attributes for writing to netcdf
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !vert. integrated 2D biomass
      character(len=10) :: b2d_his_name   = 'algae_b2d'
      character(len=12)  :: b2d_long_name  = '2D biomass'
      character(len=10)  :: b2d_units = 'g-dry/m^2'


      !3D biomass
      character(len=7) :: b_his_name   = 'algae_b'
      character(len=7)  :: b_long_name  = 'biomass'
      character(len=10)  :: b_units = 'g-dry/m^3'

      !leaf-area density
      character(len=9) :: lad_his_name = 'algae_lad'
      character(len=17)   :: lad_long_name = 'leaf-area density'
      character(len=10) :: lad_units     = '1/m'

      !algae length
      character(len=15) :: length_his_name = 'algae_length'
      character(len=17) :: length_long_name = 'algae length'
      character(len=10) :: length_units     = 'm'


      !Vertical distribution function
      character(len=11) :: shape_his_name = 'algae_shape'
      character(len=20)   :: shape_long_name = 'biomass distribution'
      character(len=10) :: shape_units     = '1/m'

      !NPP (corrector)
      character(len=11) :: npp_his_name = 'algae_npp'
      character(len=25)   :: npp_long_name = 'net primary productivity'
      character(len=15) :: npp_units     = 'g-dry / m^2 / sec'

      ! Losses ('natural')
      character(len=11) :: loss_his_name = 'algae_loss'
      character(len=25)   :: loss_long_name = 'losses'
      character(len=15) :: loss_units     = 'g-dry / m^2 / sec'

      ! harvest
      character(len=11)   :: harv_his_name = 'algae_harv'
      character(len=25)   :: harv_long_name = 'harvest'
      character(len=15)   :: harv_units     = 'g-dry / m^2 / sec'


      !Nitrate uptake
      character(len=13) :: uptno3_his_name = 'algae_uptno3'
      character(len=20)   :: uptno3_long_name = 'nitrate uptake'
      character(len=20) :: uptno3_units     = 'mmol N / m^3 / sec'

      !Exudation loss
      character(len=13) :: exuL_his_name = 'algae_exuL'
      character(len=30)   :: exuL_long_name = 'exudation loss'
      character(len=20) :: exuL_units     = 'mmol N / m^3 / sec'

      !spom loss
      character(len=13) :: spomL_his_name = 'algae_spomL'
      character(len=30)   :: spomL_long_name = 'small particulate loss'
      character(len=20) :: spomL_units     = 'mmol N / m^3 / sec'

      !frond loss
      character(len=13) :: frondL_his_name = 'algae_frondL'
      character(len=20)   :: frondL_long_name = 'large particulate loss'
      character(len=20) :: frondL_units     = 'mmol N / m^3 / sec'

      !3D PAR
      character(len=13) :: par_his_name = 'PAR'
      character(len=20)   :: par_long_name = '3D PAR'
      character(len=20) :: par_units     = 'W/m^2'

      !REMIN POC large
      character(len=20) :: reminpoclrg_his_name = 'remin_lrgpocalgae'
      character(len=20) :: reminpoclrg_long_name = 'remin of frond POC'
      character(len=20) :: reminpoclrg_units     = 'mmol C / m^3 / sec'


      !REMIN POC small
      character(len=20) :: reminpocsm_his_name = 'remin_smpocalgae'
      character(len=20) :: reminpocsm_long_name = 'remin of small POC'
      character(len=20) :: reminpocsm_units     = 'mmol C / m^3 / sec'


      !REMIN POC hrv (processed)
      character(len=20) :: reminpochrv_his_name = 'remin_hrvpocalgae'
      character(len=20) :: reminpochrv_long_name = 'remin of hrv POC'
      character(len=20) :: reminpochrv_units     = 'mmol C / m^3 / sec'


      !REMIN Kelp DOC
      character(len=20) :: remindoc_his_name = 'remin_docalgae'
      character(len=20) :: remindoc_long_name = 'remin of DOC'
      character(len=20) :: remindoc_units     = 'mmol C / m^3 / sec'


     
      !REMIN Kelp DOC refractory
      character(len=20) :: remindocref_his_name = 'remin_docrefalgae'
      character(len=30) :: remindocref_long_name = 'remin of DOCrefract'
      character(len=20) :: remindocref_units     = 'mmol C / m^3 / sec'






      ! Netcdf outputting for use in wrt_mag subroutine(s)
      ! -----------------
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_mag=0
      integer,save :: navg_mag = 0               ! number of samples in average
      integer :: record_avg = nrpf_avg           ! Triggers making of initial file
      integer :: record_his = nrpf_his
      real,save :: output_time_his = 0
      real,save :: output_time_avg = 0
     

      !printing node for debugging
      integer :: node_print=-1
  
      ! ---------------------
      ! Public functios
      ! ----------------------
      public set_mag_frc !called from set_forces.F
      public step_mag    ! called from step3d_t_iso.F (before BEC/MARBL)
      public init_arrays_mag ! called from init_arrays.F
      public wrt_mag_from_main !called from main.F

      ! =======================
      ! END OF MODULE PRE-AMBLE
      ! =======================

      contains

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


      !Initialize MAG arrays (called from init_arrays.F)
      subroutine init_arrays_mag![
      implicit none
      
      !print *, 'HELLO INIT MAG'
      !ncforce for prescribed forcing
      if (frc_biomass) allocate(nc_biomass%vdata(GLOBAL_2D_ARRAY,2) )
      if (frc_sflxpockelp) allocate(nc_sflxkelp%vdata(GLOBAL_2D_ARRAY,2) )
      
      !Read in forcing
      allocate( sflx_pockelp_mass(GLOBAL_2D_ARRAY) )
      sflx_pockelp_mass=0.


      !Standard MAG (algae) variables
      allocate( algae_age(GLOBAL_2D_ARRAY) )
      algae_age=0.
      allocate( algae_b2d(GLOBAL_2D_ARRAY) )
      algae_b2d = 0.
      allocate( algae_b(GLOBAL_2D_ARRAY,N) )
      algae_b=0.
      allocate( algae_lad(GLOBAL_2D_ARRAY,N) )
      algae_lad=0.
      allocate( algae_length(GLOBAL_2D_ARRAY) )
      algae_length=0.
      allocate( algae_length_sigcorr(GLOBAL_2D_ARRAY) )
      algae_length_sigcorr=0.
      allocate( algae_shape(GLOBAL_2D_ARRAY,N))
      algae_shape=0.
      allocate( algae_Nf(GLOBAL_2D_ARRAY,N))
      algae_Nf=0.
      allocate( algae_mask(GLOBAL_2D_ARRAY) )
      algae_mask=0.
      allocate( algae_npp(GLOBAL_2D_ARRAY) )
      algae_npp=0.
      allocate( algae_loss(GLOBAL_2D_ARRAY) )
      algae_loss=0.
      allocate( algae_harv(GLOBAL_2D_ARRAY) )
      algae_harv=0.
 
      !Converted POC kelp surface flux to BEC-compatible units
      allocate(pockelp_prod_surf(GLOBAL_2D_ARRAY))
      pockelp_prod_surf=0.


      !BEC INTERACTION VARIABLES
      allocate( algae_uptno3(GLOBAL_2D_ARRAY,N) )
      algae_uptno3 = 0.
      allocate(algae_exuL(GLOBAL_2D_ARRAY,N) )
      algae_exuL = 0.
      allocate(algae_spomL(GLOBAL_2D_ARRAY,N) )
      algae_spomL =0.
      allocate(algae_frondL(GLOBAL_2D_ARRAY,N) )
      algae_frondL =0.
      allocate(algae_harvL(GLOBAL_2D_ARRAY,N) )
      algae_harvL =0.
      allocate(PAR_3d_mag(GLOBAL_2D_ARRAY,N) )
      PAR_3d_mag =0.
      allocate(remin_lrgpocalgae(GLOBAL_2D_ARRAY,N) )
      remin_lrgpocalgae =0.    
      allocate(remin_hrvpocalgae(GLOBAL_2D_ARRAY,N) )
      remin_hrvpocalgae =0.
      allocate(remin_smpocalgae(GLOBAL_2D_ARRAY,N) )
      remin_smpocalgae =0.
      allocate(remin_docalgae(GLOBAL_2D_ARRAY,N) )
      remin_docalgae =0.
      allocate(remin_docrefalgae(GLOBAL_2D_ARRAY,N) )
      remin_docrefalgae =0.
             
      
  
      !Averages for writing
      allocate( algae_b2d_avg(GLOBAL_2D_ARRAY) )
      algae_b2d_avg=0.
      allocate( algae_b_avg(GLOBAL_2D_ARRAY,N) )
      algae_b_avg=0.
      allocate( algae_lad_avg(GLOBAL_2D_ARRAY,N) )
      algae_lad_avg=0.
      allocate( algae_length_avg(GLOBAL_2D_ARRAY) )
      algae_length_avg=0.
      allocate( algae_shape_avg(GLOBAL_2D_ARRAY,N) )
      algae_shape_avg=0.
      allocate( algae_uptno3_avg(GLOBAL_2D_ARRAY,N) )
      algae_uptno3_avg = 0.
      allocate( algae_npp_avg(GLOBAL_2D_ARRAY) )
      algae_npp_avg = 0.
      allocate( algae_loss_avg(GLOBAL_2D_ARRAY) )
      algae_loss_avg = 0.
      allocate( algae_harv_avg(GLOBAL_2D_ARRAY) )
      algae_harv_avg = 0.
      allocate( algae_exuL_avg(GLOBAL_2D_ARRAY,N) )
      algae_exuL_avg = 0.
      allocate( algae_spomL_avg(GLOBAL_2D_ARRAY,N) )
      algae_spomL_avg = 0.
      allocate( algae_frondL_avg(GLOBAL_2D_ARRAY,N) )
      algae_frondL_avg = 0.
      allocate( algae_harvL_avg(GLOBAL_2D_ARRAY,N) )
      algae_harvL_avg = 0.
      allocate( remin_lrgpocalgae_avg(GLOBAL_2D_ARRAY,N) )
      remin_lrgpocalgae_avg = 0.
      allocate( remin_hrvpocalgae_avg(GLOBAL_2D_ARRAY,N) )
      remin_hrvpocalgae_avg = 0.
      allocate( remin_smpocalgae_avg(GLOBAL_2D_ARRAY,N) )
      remin_smpocalgae_avg = 0.
      allocate( remin_docalgae_avg(GLOBAL_2D_ARRAY,N) )
      remin_docalgae_avg = 0.
      allocate( remin_docrefalgae_avg(GLOBAL_2D_ARRAY,N) )
      remin_docrefalgae_avg = 0.
      allocate( PAR_3d_mag_avg(GLOBAL_2D_ARRAY,N) )
      PAR_3d_mag_avg = 0.


      !Bottom depth/sigma level of algae arrays
      allocate(algae_kbot(GLOBAL_2D_ARRAY))
      algae_kbot=0
      allocate(algae_zbot(GLOBAL_2D_ARRAY))
      algae_zbot=0


#ifdef ALGAE_DRAG
      allocate( vegdrag_u(GLOBAL_2D_ARRAY,N) )
      allocate( vegdrag_v(GLOBAL_2D_ARRAY,N) )
      vegdrag_u=0.
      vegdrag_v=0.
#endif


     


      end subroutine init_arrays_mag!]


      !----------------------------------------
      subroutine set_mag_frc ![
      ! Read in MAG forcing
      ! This is used to read in a biomass as a function of time 
      ! and/or a surface flux of POCkelp (to simulate biomass sinking)
      
      implicit none
      
      integer :: ierr
 
      
      integer,dimension(4)   :: start
#ifdef ANA_MAG_INIT
      ! custom params for analytical mag IC
      integer :: i,j
      real    :: b2d0, age0, zbot0
                 ! farm boundaries on xr,yr grid
      real    :: x0_klp,y0_klp, x1_klp,y1_klp
#endif

      if (FIRST_TIME_STEP) then
         kelp_initialized=.false.
      endif


      !Turn on kelp after some prescribed time
      if (time.ge.twait_kelp) then
      !if (FIRST_TIME_STEP) then !--> depracated, set twait_kelp in mag.opt
         if (init_biomass .AND. .NOT. kelp_initialized) then


#ifdef ANA_MAG_INIT
       ! Analaytical MAG initial condition
       ! fill the following 3 arrays
       ! algae_b2d --> biomass [g-dry/m^2]
       ! algae_age --> age [sec]
       ! algae_zbot --> cultivation depth [m]

       ! constants
       x0_klp = 62  * 1e3
       x1_klp = 67  * 1e3
       y0_klp = 3 * 1e3
       y1_klp = 6 * 1e3 

       b2d0=100 
       age0 = 7776000.0 !3 months 
       zbot0 = 10

       do j=1,ny
          do i=1,nx
             if (xr(i,j).ge.x0_klp.and.xr(i,j).le.x1_klp) then
                if (yr(i,j).ge.y0_klp.and.yr(i,j).le.y1_klp) then
                    algae_b2d(i,j)=b2d0
                    algae_age(i,j) = age0
                    algae_zbot(i,j) = zbot0
                endif
             endif
          enddo
       enddo
         

#else
            ierr = nf90_open(mag_init_name, nf90_nowrite, ncid)
            start=1
            !print *, 'mag_init_name', mag_init_name
            !Read in initial 2D biomass [g-dry/m^2]
            call ncread(ncid,'algae_b2d_init',algae_b2d(i0:i1,j0:j1),start)
            !Read in initial age [seconds] (2D array)
            call ncread(ncid, 'algae_age_init', algae_age(i0:i1,j0:j1),start)
            !Read in algae bottom depth (2D array)
            call ncread(ncid, 'algae_zbot', algae_zbot(i0:i1,j0:j1),start)


#endif
            ! Fill some arrays that are used
            call create_algae_mask
            call fill_algae_kbot

            !------------------------------------------
            !Compute characterists for first time-step
            ! - plant length
            ! - vertical shape function
            ! - 3d biomass
            ! - leaf-area density
            ! - fixed nitrogen content
            !------------------------------------------
            call calc_mag_char_allij
            
            !set init_biomass false after initialized
            kelp_initialized=.true. 


         endif     
      endif
      
      !Read in forcing from netcdf file
      if (frc_biomass) then
          !This would be the case where you read in a biomass but dont'
          !evolve it online 
          call set_frc_data(nc_biomass, algae_b2d, 'r')
          call calc_mag_char_allij
      endif
      if (frc_sflxpockelp) call set_frc_data(nc_sflxkelp,sflx_pockelp_mass, 'r')

      !OR DEFINE FORCING ANALYTICALL IN SOME BLOCK OF CODE HERE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ----------------------------------------------------------------------
      !Convert MAG forcing to relevant quantities/units for ROMS-BEC coupling
      !-----------------------------------------------------------------------
      
      call conv_sflx_pockelp

      !print *, 'algae_lad', algae_lad     
      !print *, 'DONE WITH SET MAG FORCE'

      !Exchange variables (is this necessary??)
      call exchange_xxx(algae_b2d,sflx_pockelp_mass)
      call exchange_xxx(algae_lad,algae_b)
      call exchange_xxx(pockelp_prod_surf) 

      end subroutine set_mag_frc!]

      !----------------------------------------


      subroutine calc_mag_char_allij
      implicit none
      integer :: i,j
      !print *,'nx=', nx
      !print *, 'ny=',ny
      do j=1,ny
         do i=1,nx
            if (algae_mask(i,j).eq.1.D0) then
               !print *, 'calc mag char ', i, j
               call calc_mag_char(i,j)
               if (node_print.eq.-1) then
                  node_print=mynode
               endif

            endif
         enddo
      enddo
      end subroutine calc_mag_char_allij
      


      subroutine fill_algae_kbot 
      implicit none
      integer :: i,j, kbot_temp
      do j=1,ny
         do i=1,nx
            call find_k(i,j,-algae_zbot(i,j),kbot_temp)
            algae_kbot(i,j) = kbot_temp
         enddo
      enddo

      end subroutine fill_algae_kbot 


      subroutine create_algae_mask
      implicit none
      !Create mask based on 2D biomass 
      where(algae_b2d.gt.0.D0)
          algae_mask=1.D0
      end where
      end subroutine create_algae_mask

      subroutine conv_sflx_pockelp
      implicit none
      ! inputs
      integer :: i,j,k
      real :: conv
      ! presumably can just work on interior points, since no algae
      ! variables will be placed on stencil that requires buffers
      ! (...but maybe for drag???!?)
      do j=1,ny
        do i=1,nx
 
            ! --------------------------------------
            !         Convert pockelp surface flux to 
            !         units of [mmol C / m^3 / sec]
            !         from     [g C / m^2 / day]
            ! --------------------------------------
            !Convert g C / m^2 / day to mmol C / m^2 / sec
                   ! g C --> mmol C   ! 1/day --> 1/sec
            conv = (1000. / 12.01) * (1./86400.)
            !divide by Hz at surface to make it a source term at surface and
            !give units of per m^3
            !print *, 'Hz(i,j,N) = ', Hz(i,j,N)
            pockelp_prod_surf(i,j) = sflx_pockelp_mass(i,j)*conv / Hz(i,j,N)

        enddo !< -- i
      enddo !<-- j
      end subroutine conv_sflx_pockelp
      


      !===========================================

      ! PROGONSTIC MAG SOLVER
      ! SOLVES dB_2D / dt = growth - losses
      ! This is called right before BEC in step3d_t_iso.F
      !==============================================

      subroutine calc_mag_char(i,j)
      !Compute algae characteristics
      !   - plant length
      !   - vertical shape function
      !   - 
      
      implicit none

      integer,intent(in) :: i,j
      integer :: k

      !----------------------
      !Calculate algae_length (based on z-integrated biomass)
      !-----------------------
      call calc_algaelength(i,j)
      !print *, 'algae_length(i,j) in mag_char', algae_length(i,j)
      ! -----------------------------------------------------
      ! Convert vert. integral of biomass to a vert. profile
      ! ----------------------------------------------------

      ! computes algae_shape(i,j)
      call calc_shape(i,j)
      !Now distribute biomass
      do k=1,N
         algae_b(i,j,k) = algae_b2d(i,j) * algae_shape(i,j,k)
      enddo

      ! ---------------------------------------
      ! Convert biomass(z) to leaf-area density
      ! ---------------------------------------
      call calc_LAD(i,j)

      !Compute fixed nitrogen for BEC shading
      algae_Nf(i,j,:) = Qmin * algae_b(i,j,:)

      end subroutine calc_mag_char



      !----------------------------------------

      !------------------------------------------
      subroutine step_mag
      ! time-step 2D biomass

      implicit none

      integer :: i,j,k
      real(8)    :: b2d_new,loss,npp_pred,npp_corr, harv
      real    :: no3_avg, PAR_avg, corr_thick
      real(8)    :: uptno3_vint, loss_frond, loss_spom, loss_exu
      real    :: MM_no3(N)
      integer :: k_top
      logical :: print_stuff=.false.
      !integer :: node_print=-1
#ifdef MAG_EVOLVE
      if(mynode==0) then
        write(*,'(7x,A)') 'step_mag:: evolving macroalage with MoMAG'
      endif

      do i=1,nx
         do j=1,ny
            !Only do where there is kelp
            if (algae_mask(i,j).eq.1) then

               !----------------------------------
               ! Compute vert. averages of no3,PAR
               ! over the length of algae
               !---------------------------------
               !     Vertically integrate from cult. depth to top of
               !     plant
               !  note to self: checked with print that tracer array
               !  here is in actual units of tracer (and not normalized
               !  by dz), so we are fine calling it like this
               call vert_int_algae(i,j,t(i,j,:,nstp,iNO3),no3_avg) !make sure of nstp
               call vert_int_algae(i,j,PAR_3d_mag(i,j,:),PAR_avg) 
               !Divide by (rho-level respectful) length to get vertical average
               if (print_stuff) then
                 if (mynode.eq.node_print) then
                   if (mod(algae_age(i,j),43200.0).eq.0) then
                         call find_k(i,j,z_w(i,j,algae_kbot(i,j)-1) + algae_length(i,j),k_top)
                          !call find_k(i,j,-algae_zbot(i,j)+algae_length(i,j),k_top)
                          print *, '------------------------------'
                          print *, ' '
                          print*, 'no3_int=', no3_avg
                          print *, 'no3(kbot)=', t(i,j,algae_kbot(i,j),nstp,iNO3)
                          print *, 'no3(ktop)=', t(i,j,k_top,nstp,iNO3)
                          print *, 'Hz(kbot)=', Hz(i,j,algae_kbot(i,j))
                          print *, 'Hz(ktop)=', Hz(i,j,k_top)
                          print *, 'kbot=', algae_kbot(i,j)
                          print *,'ktop=', k_top
                          print *, 'z_r(kbot)=', z_r(i,j,algae_kbot(i,j))
                          print *, 'z_r(ktop)=', z_r(i,j,k_top)
                          print *, 'z_w(kbot)=', z_w(i,j,algae_kbot(i,j))
                          print *, 'z_w(ktop)=', z_w(i,j,k_top)
                          print *, 'z_r length=', z_r(i,j,k_top) - z_r(i,j,algae_kbot(i,j))
                          print *, 'algae_length=', algae_length(i,j)
                          print*, 'algae_length_sigcorr=',algae_length_sigcorr(i,j)
                          print *, 'Hz(kbot)=', Hz(i,j,algae_kbot(i,j))
                          print *, ' '
                          print *, '------------------------------'
                         !print_stuff=.false.
                    endif
                  endif
               endif 
               no3_avg = no3_avg / algae_length_sigcorr(i,j)
               PAR_avg = PAR_avg / algae_length_sigcorr(i,j) 

#ifdef CONST_MAG_FRC
               !if (mynode.eq.0) then
                  ! check tracer vaule (multiply by dz or not)
               !   print *, 'no3(1,1,N)=', t(i,j,N,nstp,iNO3)
               !endif
               !TESTING MoMAG with const. forcing
               no3_avg = 3.0 !mmol N / m^3
               PAR_avg = 50.0  !W/m^2
#endif

               !Check for negative values
               if (no3_avg.lt.0) then
                  no3_avg=0.0
               endif        
               if (PAR_avg.lt.0) then
                  PAR_avg=0.0
               endif       
               if (print_stuff) then
                  if (mynode.eq.node_print) then
                     if (mod(algae_age(i,j),43200.0).eq.0) then
                        print *, '======================'
                        print *, ' '
                        print *, 'no3_avg=', no3_avg
                        print *, 'PAR_avg=', PAR_avg
                        print *, ' '
                        print *, '======================'
                      endif
                  endif
               endif
               !----------------------------------
               ! Compute npp (net primary productivity) term
               !    
               !    - first compute the `predictor' npp based on MoMAG
               !      rules
               !   - then correct npp based on the vertical integral of
               !     uptake
               !----------------------------------
               
               !----------------------------------
               !          PREDICTOR NPP
               !-----------------------------------
               if (algae_age(i,j).lt.juve_age) then
                  !print *, 'JUVE NPP'
                  ! juvenile is just prescribed growth rate
                  npp_pred = g0_juve * algae_b2d(i,j)
               else
                  !print *, 'MATURE NPP'
                  !mature cares about no3, par and age
                  npp_pred = algae_b2d(i,j) * ( NPP_coeff_no3 * no3_avg 
     &                                     + NPP_coeff_par * PAR_avg
     &                                     + NPP_coeff_age *(algae_age(i,j)-juve_age) )
              endif               
              
              !Don't let npp_pred exceed cap (should only be relevant
              !for extreme situations)
              npp_pred = min(npp_pred,NPP_max)

#ifdef DPD_REGRESSION
               npp_pred = algae_b2d(i,j) * ( dpd_coeff_intercept  
     &                                        +  dpd_coeff_no3 * no3_avg
     &                                        +  dpd_coeff_PAR * PAR_avg)
#endif


#ifdef FAKE_NPP_LOSS
             !WHILE WE WORK OUT COEFFICIENT UNITS JUST HAVE A HARDCODED
             !NPP FOR TESTING
             npp_pred = g0_fake * algae_b2d(i,j)         
#endif

#ifdef MAG_GL_AMP
               !Artificially amplify growth based on coefficient in
               !mag.opt (for ensemble of growth scenarios)
               npp_pred = npp_pred* alpha_growth
#endif
              !----------------------------------
               !          CORRECTOR NPP
               !        compute vert. variable no3 uptake
               !        then integrate back to get NPP
               !-----------------------------------
               
               do k=1,N
                   !Micahelis Menten uptake coefficient -->actually
                   !looks ot result in  big damping of npp
                   !MM_no3(k) = t(i,j,k,nstp,iNO3) / ( t(i,j,k,nstp,iNO3)+ Kno3_mag)
                   
                   !Heaviside function (trying it out)
                   if (t(i,j,k,nstp,iNO3).gt.0) then
                      MM_no3(k)=1.
                   else
                      MM_no3(k)=0.
                   endif

               enddo

               !This is the uptake that BEC will see
               algae_uptno3(i,j,:) = MM_no3(:) * algae_shape(i,j,:) *(Qmin / MW_N) * npp_pred
              

               !Now, correct NPP based on the vertical integral of the
               !uptake over the kelp thickness
               call vert_int_algae(i,j,algae_uptno3(i,j,:),uptno3_vint)
               !call vert_int_kbot_to_surf(i,j,algae_uptno3(i,j,:),uptno3_vint)
               
               !corrected NPP is then
               npp_corr = (MW_N / Qmin) * uptno3_vint
#ifdef NO_NPP_CORR
               npp_corr = npp_pred
#endif
               ! assign for writing
               algae_npp(i,j) = npp_corr
               

               !print *, 'NPP_CORR - NPP_PRED=', npp_corr - npp_pred
               !------------------------------------
               ! Compute MoMAG loss term
               !-------------------------------------
               if (algae_age(i,j).lt.juve_age) then
                  !NO LOSSES FOR 'BABY' KELP
                  loss = 0.D0
               else
                  !MoMAG regression model
                  !print *, 'par loss = ', algae_b2d(i,j) * L_coeff_par * PAR_avg
                  !print *, 'age loss = ', algae_b2d(i,j) *  L_coeff_age * (algae_age(i,j) - juve_age)
                  !print *, 'total = ', algae_b2d(i,j) * L_coeff_par * PAR_avg 
! &                    +  algae_b2d(i,j) *  L_coeff_age *  (algae_age(i,j) - juve_age)
                  !loss = algae_b2d(i,j) * (L_coeff_par * PAR_avg + L_coeff_age * (algae_age(i,j) - juve_age) )

                  ! New MoMAG has no3, par, and age dependent loss                  
                  loss = algae_b2d(i,j) * ( L_coeff_no3 * no3_avg 
     &                                     + L_coeff_par * PAR_avg
     &                                     + L_coeff_age *(algae_age(i,j)-juve_age) )
 
               endif

               ! assign for writing 
               algae_loss(i,j) = loss



               !Harvest
               !        biomass based threshold 
               !        user can implement whatever scenario they want
               !if (algae_b2d(i,j).gt.harvest_thresh_b2d) then
               if (algae_age(i,j).ge.harvest_t0.and.algae_age(i,j).lt.harvest_t1) then
                  harv = algae_b2d(i,j) * harvest_rate
               else
                  harv=0.D0
               endif
               algae_harv(i,j) = harv 
 
               !-------------------------------
               !  TIME STEP w/ forward euler 
               !---------------------------------------
               b2d_new = algae_b2d(i,j) + dt * (npp_corr -  loss - harv)
               ! leave like this for now for some checks
               if (print_stuff) then
                if (mynode.eq.node_print) then
                  if (mod(algae_age(i,j),43200.0).eq.0) then
                     print *, '***************************'
                     print *, ' '
                     print *, '    npp = ', npp_corr
                     print*,  '    npp_pred=', npp_pred
                     print*,  '    pred - corr ', npp_pred - npp_corr
                     print *, '    loss = ', loss
                     print *, '    db = ', npp_corr + loss
                     print *, ' '
                     print *, '***************************'
                  endif
                endif
               endif
               !---------------------------------------
               ! PARTITION/ASSIGN BEC COUPLING
               ! - before updating length, shape, etc.               
               !--------------------------------------
               !arbitrarily divide losses equally into 1/3s
               
               ! frond, spom, exu are NATURAL losses (harv separate)

               ! LOSS IS POSITIVE, NO MORE FLIPPING SIGNS
               !large frond (natural losses)
               loss_frond = loss * coeff_frondL 
               !small particulates
               loss_spom  = loss * coeff_spomL
               !exudation of DOM
               loss_exu   = loss * coeff_exuL

               !Now convert to 3D coupling terms for BEC
               ! note that loss_ already has biomass units (see code
               ! above)
               ! this is different than notation in paper where L is a
               ! rate
               algae_exuL(i,j,:) = (Qmin/MW_N)* algae_shape(i,j,:) * loss_exu
               algae_frondL(i,j,:) = (Qmin/MW_N)* algae_shape(i,j,:) * loss_frond
               algae_spomL(i,j,:) = (Qmin/MW_N)* algae_shape(i,j,:) *loss_spom

               !a bit odd to scale the harvest by the shape of the farm,
               ! but unsure of a better/self-consistent alternative 
               algae_harvL(i,j,:) = (Qmin/MW_N)* algae_shape(i,j,:) * harv



               !------------------------------
               ! Update kelp characteristics with new biomass
               !  - compute plant length
               !  - compute shape function
               !  - vert. distribute into 3d biomass
               !  - compute LAD
               !  - compute Nf for BEC shading
               !-------------------------------
               !update b2d array
               !                   restrict negative biomass (should
               !                   never happen...)
               algae_b2d(i,j) = max(0.D0,b2d_new)

               call calc_mag_char(i,j) 
               
               !Update age
               algae_age(i,j) = algae_age(i,j) + dt
               !print *, 'age=', algae_age(i,j), ' sec'
               !print *, 'age=', algae_age(i,j)/86400., ' days'
          endif !<---algae mask
        enddo
      enddo
      
      !Exchange LAD for translation to u,v points in drag
      call exchange_xxx(algae_lad)
      !call exchange_xxx(PAR_3d_mag)
      !call exchange_xxx(algae_uptno3)
      !call exchange_xxx(algae_exuL,algae_frondL,algae_spomL)
      !call exchange_xxx(remin_pocalgae)
#endif 
      end subroutine step_mag


      !======================================================
      ! 
      !                   Algae characteristic functions
      !
      !======================================================

      subroutine calc_shape(i,j)
      
      !------------------------------------
      !Calculate vert distribution function
      !------------------------------------
      
      implicit none
      integer,intent(in):: i,j
      integer :: k,k0,k_top
      real :: ztop, fexc,fcan, Bm_noc
      real :: zloc(N),zz(N),zmask(N), Bm_can(N)
      real :: temp_shape(N)
      real :: temp_shape_int


      !-----------------------------------------------------------
      !Implementation of D. Bianchi 'constant plus' shape function
      !-----------------------------------------------------------

      !top of water column starting from algae_zbot (which is positive)
      ztop = zeta(i,j,nstp) + algae_zbot(i,j) 
      !ztop = 0 - algae_zbot(i,j)
      !Plant coordinate systems in z
      do k=1,N
         !define a local depth coordinate
         zloc(k) = z_r(i,j,k) + algae_zbot(i,j)
         !normalized depth coordinate
         zz(k) = zloc(k) / algae_length(i,j)
         !plant mask in vertical (1=plant,0=water)
         zmask(k) = 1.0
         if (zz(k)<0) then
            zmask(k) = 0.0
         endif
         if (zz(k)>1) then
            zmask(k)=0.0
         endif
      enddo
      
      ! Define fraction of plant excess (length - ztop) 
      ! relative to maximum possible plant excess (Hmax - ztop)
     !print *, 'Hmax=', Hmax
      if (algae_length(i,j).le.ztop) then
         fexc = 0.0
      else
          fexc = (algae_length(i,j) - ztop) / (Hmax - ztop)
      endif
 
      fcan  = cmax_shape * fexc
      !No canopy baseline (constant)
      Bm_noc = 1 / min(algae_length(i,j), ztop)

      do k=1,N
        !Canopy forming section
        Bm_can(k) = 1 / (Hc * ( 1- exp(-ztop/Hc)))*exp((zloc(k)-ztop)/Hc) 
        temp_shape(k) = (1 - fcan) * Bm_noc + fcan * Bm_can(k)
        !mask
        temp_shape(k) = temp_shape(k) * zmask(k)
      enddo
      !Get sigma level of top of plant
      !k_top = algae_ktop(i,j)
      call find_k(i,j,z_w(i,j,algae_kbot(i,j)-1) + algae_length(i,j),k_top)
      !call find_k(i,j,-algae_zbot(i,j) + algae_length(i,j),k_top)
     

      if (const_Bm_shape) then
         !Make a dummy function that is constant in k
         temp_shape(:)=0.
         do k=k0,k_top
              temp_shape(k)=1.
         enddo
      endif
      !print *, 'shape(k0+1)', temp_shape(k0+1)

      !Integrate from cultivation depth to top of plant     
      call vert_int_algae(i,j,temp_shape,temp_shape_int)
       
      !Normalize
      !         This is not formally necessary, as the analytical
      !         function has unit integral, but is enforced in case of
      !         numerical errors 
      algae_shape(i,j,:) = temp_shape(:) / temp_shape_int
      end subroutine calc_shape

      subroutine calc_algaelength(i,j)
      !==============================================
      ! Calculate kelp length based on vertical integral of biomass
      !==============================================
      implicit none
      !inputs
      integer, intent(in):: i,j  
      real :: temp_length
      integer :: ktop
      !Hmax, Kh set in mag_param.opt
      temp_length = (Hmax * algae_b2d(i,j) ) 
     &                   / (Kh + algae_b2d(i,j) )
      if (satellite_frc) then         
         !If forcing with remotely sensed data, kelp reaches surface by
         !definition
         algae_length(i,j) = max(algae_zbot(i,j),temp_length)
         !print *, 'KELPWATCH FRC, temp_length = ', temp_length, '  zbot  = ', algae_zbot(i,j)
      else
         algae_length(i,j) = temp_length
      endif


      !---------------------------------------------------------
      !           NOTES on vertical discretization of kelp
      !   We place kelp (algae_b) at vertical rho-levels, to make the
      !   coupling to ROMS-BEC straightfoward, as all the variables the
      !   kelp will interact with are at vertical rho-levels (..except
      !   Kv).  So, in the limitng case that the kelp length is less
      !   than one grid-cell thick (kbot=ktop), the kelp is assumed to
      !   fill a vertical grid cell from z_w(i,j,kbot-1) to z_w(i,j,kbot) with
      !   a thickness Hz(kbot). 
      !   In the more general case, we assume the kelp extends from the
      !   w-level interace below algae_zbot to the w-level interface
      !   above ktop, which is just z_w(i,j,ktop) - z_w(i,j,kbot); that
      !   is, when we compute a vertical integral of something the kelp
      !   'sees' it is over this thickness
      ! 
      !   Additionally, the vertical thickness of the kelp cannot be
      !   greater than the water-column thickness (up from the
      !   cultivation depth); in the case of a canopy, the kelp length
      !   can be longer, but the thickness cannot (ktop cannot be
      !   greater than N). So, we do a check on ktop below (but find_k
      !   should not give anything greater than N anyways)
      !
      !   These choices and implementation of the kelp on the vertical
      !   grid DO imply that solutions are sensitive to vertical
      !   resolution and you probably want small dz (high res) where you
      !   have kelp (...so stretch there)
      !---------------------------------------------------------
     
      !Find sigma level corresponding to top of kelp

      call find_k(i,j,z_w(i,j,algae_kbot(i,j)-1) + algae_length(i,j),ktop)
 
      !call find_k(i,j,-algae_zbot(i,j) + algae_length(i,j),ktop)
      !Ensure ktop is not greater than N (shouldn't happen, but its a
      !check for canopy-forming)

      ktop = min(ktop,N)
      !Compute `corrected' length that corresponds to kelp
      !discretization on sigma levels
      
      ! top of platn is at w-level k-top, bottom is a at w-level kbot-1
      ! b/c we assume the kelp at a rho-point fills the entire grid cell
      ! (from w-level to w-level interface)

      !this the same as summing Hzs from kbot to ktop
      algae_length_sigcorr(i,j) = z_w(i,j,ktop) - z_w(i,j,algae_kbot(i,j)-1)
      !algae_length_sigcorr(i,j) = z_r(i,j,ktop) - z_r(i,j,algae_kbot(i,j))



      !if (ktop.eq.algae_kbot(i,j)) then
      !   algae_length_sigcorr(i,j) = Hz(i,j,ktop) 
      !else
      !  algae_length_sigcorr(i,j) = z_r(i,j,ktop) - z_r(i,j,algae_kbot(i,j))
      !endif

      !print *, 'length - sig corr length=', algae_length(i,j) - algae_length_sigcorr(i,j) 
      end subroutine calc_algaelength

      subroutine calc_LAD(i,j)
      !==============================================
      ! Calculate leaf-area density 
      !==============================================
      implicit none
      !inputs
      integer, intent(in):: i,j  

      ! Convert biomass to leaf-area density
      ! lad = b / 0.094 * sa2b0 ,where sa2b0
      ! is a surface area to biomass ratio (defined in mag_param.opt)
      !         in the old version, sa2b was defined in a piecewise way,
      !         but this would result in an `ugly' (=discontinous)
      !         vertical structure of LAD, regardless of the vertical
      !         structure of algae_b (biomass)
      !      
      !         with picking one sa2b0=something, you let algae_lad take
      !         the vertical structure of algae_b (from algae_shape
      !         [1/m] vertical distribution function)

  
      !                     [g-dry/m^3] /  [g-wet/g-dry] * [m^2 / g-wet]
      !                     = 1 / m
      algae_lad(i,j,:) = (algae_b(i,j,:) / dry_wet) * sa2b0 

      end subroutine calc_LAD

      !======================================================
      ! 
      !                   UTILITY FUNCTIONS
      !
      !======================================================
      subroutine find_k(i,j,z_look,k_out)
      !================================
      !Find sigma level corresponding
      ! to a specified depth (z_look)
      !================================
      implicit none
      real, intent(in) :: z_look !should be negative
      integer, intent(in) :: i,j
      integer, intent(out) :: k_out

      integer :: k
      real    :: min_diff, diff
      !real :: diff(N)
      !diff = abs(z_r(i,j,:) - z_look)
      !k_out = minloc(diff,1)
     
       !Courtesy of chatgpt 
       min_diff=1.0e20
       k_out=1 !default
       do k=1,N
          diff = abs(z_r(i,j,k) - z_look)   ! Compute the absolute difference
          if (diff < min_diff) then
             min_diff = diff             ! Update the minimum difference
             k_out = k                   ! Update the closest level index
          end if
       enddo
      end subroutine find_k


      subroutine vert_int_kbot_to_surf(i,j,var_in,var_out)
      ! Compute vertical integral from cultivation
      ! depth to surface
      implicit none
      integer, intent(in) :: i,j
      real, dimension(:),intent(in) :: var_in(N)
      real, intent(out) :: var_out
      integer :: k,k0
      !k0=algae_kbot(i,j)  !kbot!(i,j) 
      var_out = 0
      !call find_k(i,j,algae_zbot*-1,k0)
      do k=algae_kbot(i,j),N
         var_out = var_out  
     &            + var_in(k) * Hz(i,j,j)
      enddo
      !do k=algae_kbot(i,j),N-1
      !   var_out=var_out
     !&          + 0.5 * (var_in(k+1) + var_in(k))
     !&          * (z_w(i,j,k+1) -z_w(i,j,k))
     ! enddo

      end subroutine vert_int_kbot_to_surf 
      
      subroutine vert_int_algae(i,j,var_in,var_out)
      ! Compute vertical integral from cultivation
      ! depth to top of plant
      implicit none
      integer, intent(in) :: i,j
      real, dimension(:),intent(in) :: var_in(N)
      real, intent(out) :: var_out !,algae_corr_thickness
      integer :: k,k0,ktop
     
      var_out = 0
      !find bottom (k=k0) 
      !call find_k(i,j,algae_zbot(i,j)*-1,k0)
      !fidn top of plant (k=ktop)
      
      !use w-level kbot to get top, assuming kelp fills grid cell
      call find_k(i,j,z_w(i,j,algae_kbot(i,j)-1) + algae_length(i,j),ktop)
      !call find_k(i,j,-algae_zbot(i,j) + algae_length(i,j),ktop)

      do k=algae_kbot(i,j),ktop
         var_out = var_out 
     &           + var_in(k) * Hz(i,j,k) 
      enddo


      !if (algae_kbot(i,j).eq.ktop) then
      !   var_out = var_in(ktop) * Hz(i,j,ktop)
      !else
      !   do k=algae_kbot(i,j),ktop
      !      var_out=var_out
     !&          + 0.5 * (var_in(k) + var_in(k+1))
     !&          * (z_r(i,j,k+1) -z_r(i,j,k))
     !    enddo
     ! endif
      !output 'corrected' length, reflecting vert. coord system
     
      !algae_corr_thickness = z_r(i,j,ktop) - z_r(i,j,algae_kbot(i,j))
      !print *, 'length - corr thickness =' algae_length(i,j) - algae_corr_thickness
      end subroutine vert_int_algae
      !======================================================





      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                 WRITING OUTPUT SUBROUTINES
      !                 these follow functionality in
      !                 surf_flux.F

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! NEW I/O Sept 2024
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ----------------------------------------------------------------------
      subroutine wrt_mag_from_main(init)  ![
      implicit none
! import/export
      logical,optional       :: init
       
      call wrt_mag
      ! initialization writting
      !if (wrt_file_his.and.present(init)) then
      !   call wrt_mag(.false.) ! write his data
      !else
      ! local
      !   if (wrt_file_his)    call wrt_mag(.false.) ! write his data
      !   if (wrt_file_avg)    call wrt_mag(.true.) ! write avg data
      !endif

      end subroutine wrt_mag_from_main  !]

      !---------------------------------------------------------------------------
      subroutine wrt_mag  ![
      ! write surface flux variables to netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      character(len=99),save :: fname_his,fname_avg
      integer                :: ierr
      logical,save           :: first_step_avg=.true.
      logical,save           :: first_step_his=.true. 




      ! WRITE AVERAGE
      if (wrt_file_avg) then
           call calc_avg_mag                          ! don't include t=0 in averaging
           output_time_avg = output_time_avg + dt     ! only start count after first timestep
        if (output_time_avg>=output_period_avg) then  ! time for an output
          if (mod(record_avg,nrpf_avg)==0) then
             call create_mag_file(fname_avg,.true.)
             record_avg = 0
          endif
          record_avg = record_avg + 1
          output_time_avg = 0
          navg_mag = 0
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_avg/))
          !call ncwrite(ncid,'ocean_time',(/t_avg_mag/),(/record_avg/))
          if (mynode==0) then
             write(*,'(7x,A,1x,A)')
     &       'wrt_mag :: average'
          endif
          start=1; start(3)=record_avg
          !2D outputs
          if (wrt_avg_b2d) call ncwrite(ncid,b2d_his_name,algae_b2d_avg(i0:i1,j0:j1),start)
          if (wrt_avg_npp) call ncwrite(ncid,npp_his_name,algae_npp_avg(i0:i1,j0:j1),start)
          if (wrt_avg_loss) call ncwrite(ncid,loss_his_name,algae_loss_avg(i0:i1,j0:j1),start)
          if (wrt_avg_length) call ncwrite(ncid,length_his_name,algae_length_avg(i0:i1,j0:j1),start)
          
          !3D outputs
          start=1; start(4)=record_avg
          if (wrt_avg_b) call ncwrite(ncid,b_his_name,algae_b_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_shape) call ncwrite(ncid,shape_his_name,algae_shape_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_lad) call ncwrite(ncid,lad_his_name,algae_lad_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_uptno3) call ncwrite(ncid,uptno3_his_name,algae_uptno3_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_exuL) call ncwrite(ncid,exuL_his_name,algae_exuL_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_spomL) call ncwrite(ncid,spomL_his_name,algae_spomL_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_frondL) call ncwrite(ncid,frondL_his_name,algae_frondL_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_par) call ncwrite(ncid,par_his_name,PAR_3d_mag_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_reminpoclrg) call ncwrite(ncid,reminpoclrg_his_name,remin_lrgpocalgae_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_reminpocsm) call ncwrite(ncid,reminpocsm_his_name,remin_smpocalgae_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_remindoc) call ncwrite(ncid,remindoc_his_name,remin_docalgae_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_remindocref) call ncwrite(ncid,remindocref_his_name,remin_docrefalgae_avg(i0:i1,j0:j1,:),start)
 
          if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'mag :: wrote average, tdays =', tdays,
     &       'step =', iic, 'rec =', record_avg
          endif
      endif !-->output_time avg if
      endif !-->wrt_file_Avg
      
      ! WRITE HISTORY
      if (wrt_file_his) then
        if (.not. first_step_his) then
           output_time_his = output_time_his + dt
        endif
        first_step_his=.false.
        if (output_time_his>=output_period_his .or.
     &      output_time_his==0                 ) then  ! time for an output
          if (mod(record_his,nrpf_his)==0) then
             call create_mag_file(fname_his,.false.)
             record_his = 0
          endif
          record_his = record_his + 1
          output_time_his=0
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_his/))
          start=1; start(3)=record_his
          !2D outputs
          if (wrt_b2d) call ncwrite(ncid,b2d_his_name,algae_b2d(i0:i1,j0:j1),start)
          if (wrt_npp) call ncwrite(ncid,npp_his_name,algae_npp(i0:i1,j0:j1),start)
          if (wrt_loss) call ncwrite(ncid,loss_his_name,algae_loss(i0:i1,j0:j1),start)
          if (wrt_harv) call ncwrite(ncid,loss_his_name,algae_harv(i0:i1,j0:j1),start)
          if (wrt_length) call ncwrite(ncid,length_his_name,algae_length(i0:i1,j0:j1),start)
 
          !3D outputs
          start=1; start(4)=record_his
          if (wrt_b) call ncwrite(ncid,b_his_name,algae_b(i0:i1,j0:j1,:),start)
          if (wrt_shape) call ncwrite(ncid,shape_his_name,algae_shape(i0:i1,j0:j1,:),start)
          if (wrt_lad) call ncwrite(ncid,lad_his_name,algae_lad(i0:i1,j0:j1,:),start)
          if (wrt_uptno3) call ncwrite(ncid,uptno3_his_name,algae_uptno3(i0:i1,j0:j1,:),start)
          if (wrt_exuL) call ncwrite(ncid,exuL_his_name,algae_exuL(i0:i1,j0:j1,:),start)
          if (wrt_spomL) call ncwrite(ncid,spomL_his_name,algae_spomL(i0:i1,j0:j1,:),start)
          if (wrt_frondL) call ncwrite(ncid,frondL_his_name,algae_frondL(i0:i1,j0:j1,:),start)
          if (wrt_par) call ncwrite(ncid,par_his_name,PAR_3d_mag(i0:i1,j0:j1,:),start)
          if (wrt_reminpoclrg) call ncwrite(ncid,reminpoclrg_his_name,remin_lrgpocalgae(i0:i1,j0:j1,:),start)
          if (wrt_reminpocsm) call ncwrite(ncid,reminpocsm_his_name,remin_smpocalgae(i0:i1,j0:j1,:),start)
          if (wrt_remindoc) call ncwrite(ncid,remindoc_his_name,remin_docalgae(i0:i1,j0:j1,:),start)
          if (wrt_remindocref) call ncwrite(ncid,remindocref_his_name,remin_docrefalgae(i0:i1,j0:j1,:),start)
 
          if (mynode == 0) then
             write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'mag :: wrote history, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record_his
          endif
         endif
      endif !-->wrt_file_his        
      end subroutine wrt_mag !]
      !---------------------------------------------------------------------------



      ! ----------------------------------------------------------------------
      subroutine create_mag_file(fname,avg)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ncid,ierr,varid

      if (avg) then
        call create_file('_mag_avg',fname)
      else
        call create_file('_mag',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call def_vars_mag( ncid,avg )

      ierr = nf90_enddef(ncid)

      end subroutine create_mag_file !]
      ! ----------------------------------------------------------------------


      ! ----------------------------------------------------------------------
      subroutine def_vars_mag( ncid, avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for actual or averaged variables
      implicit none

      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid, itrc
      character(len=64) :: text_lname


      !2D Biomass
      if ( ( wrt_b2d .and. .not. avg) .or. (wrt_avg_b2d.and. avg ) ) then
        if (.not. avg) text_lname=b2d_long_name
        if (      avg) text_lname='averaged '// b2d_long_name
        varid =
     &  nccreate(ncid,b2d_his_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',b2d_units)
      endif

      !3D Biomass
      if ( ( wrt_b .and. .not. avg) .or. (wrt_avg_b.and. avg ) ) then
        if (.not. avg) text_lname=b_long_name
        if (      avg) text_lname='averaged '// b_long_name
        varid =
     &  nccreate(ncid,b_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',b_units)
      endif

      !Leaf-area density
      if ( ( wrt_lad .and. .not. avg) .or. (wrt_avg_lad .and. avg ) ) then
        if (.not. avg) text_lname=lad_long_name
        if (      avg) text_lname='averaged '// lad_long_name
        varid =
     &  nccreate(ncid,lad_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',lad_units)
      endif

      !length
      if ( ( wrt_length .and. .not. avg) .or. (wrt_avg_length.and. avg ) ) then
        if (.not. avg) text_lname=length_long_name
        if (      avg) text_lname='averaged '// length_long_name
        varid =
     &  nccreate(ncid,length_his_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',length_units)
      endif


      !Biomass shape function
      if ( ( wrt_shape .and. .not. avg) .or. (wrt_avg_shape .and. avg )) then
        if (.not. avg) text_lname=shape_long_name
        if (      avg) text_lname='averaged '// shape_long_name
        varid =
     &  nccreate(ncid,shape_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',shape_units)
      endif


      ! npp
      if ( ( wrt_npp .and. .not. avg) .or. (wrt_avg_npp.and. avg ) ) then
        if (.not. avg) text_lname=npp_long_name
        if (      avg) text_lname='averaged '// npp_long_name
        varid =
     &  nccreate(ncid,npp_his_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',npp_units)
      endif

      !loss 
      if ( ( wrt_loss .and. .not. avg) .or. (wrt_avg_loss.and. avg ) ) then
        if (.not. avg) text_lname=loss_long_name
        if (      avg) text_lname='averaged '// loss_long_name
        varid =
     &  nccreate(ncid,loss_his_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',npp_units)
      endif

      !harvest 
      if ( ( wrt_harv .and. .not. avg) .or. (wrt_avg_harv.and. avg ) ) then
        if (.not. avg) text_lname=harv_long_name
        if (      avg) text_lname='averaged '// harv_long_name
        varid =
     &  nccreate(ncid,harv_his_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',npp_units)
      endif


      !Nitrate Uptake
      if ( ( wrt_uptno3 .and. .not. avg) .or. (wrt_avg_uptno3 .and. avg) ) then
        if (.not. avg) text_lname=uptno3_long_name
        if (      avg) text_lname='averaged '// uptno3_long_name
        varid =
     &  nccreate(ncid,uptno3_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',uptno3_units)
      endif

      !exudation loss
      if ( ( wrt_exuL .and. .not. avg) .or. (wrt_avg_exuL .and. avg) ) then
        if (.not. avg) text_lname=exuL_long_name
        if (      avg) text_lname='averaged '// exuL_long_name
        varid =
     &  nccreate(ncid,exuL_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',exuL_units)
      endif

      !small particulate loss
      if ( ( wrt_spomL .and. .not. avg) .or. (wrt_avg_spomL .and. avg) ) then
        if (.not. avg) text_lname=spomL_long_name
        if (      avg) text_lname='averaged '// spomL_long_name
        varid =
     &  nccreate(ncid,spomL_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',spomL_units)
      endif


      !frond particulate loss
      if ( ( wrt_frondL .and. .not. avg) .or. (wrt_avg_frondL .and. avg) ) then
        if (.not. avg) text_lname=frondL_long_name
        if (      avg) text_lname='averaged '// frondL_long_name
        varid =
     &  nccreate(ncid,frondL_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',frondL_units)
      endif



      !PAR 3d
      if ( ( wrt_par .and. .not. avg) .or. (wrt_avg_par .and. avg) ) then
        if (.not. avg) text_lname=par_long_name
        if (      avg) text_lname='averaged '// par_long_name
        varid =
     &  nccreate(ncid,par_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',par_units)
      endif



      !remin of frond POC
      if ( ( wrt_reminpoclrg .and. .not. avg) .or. (wrt_avg_reminpoclrg .and. avg) ) then
        if (.not. avg) text_lname=reminpoclrg_long_name
        if (      avg) text_lname='averaged '// reminpoclrg_long_name
        varid =
     &  nccreate(ncid,reminpoclrg_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',reminpoclrg_units)
      endif


      !remin of small POC
      if ( ( wrt_reminpocsm .and. .not. avg) .or. (wrt_avg_reminpocsm .and. avg) ) then
        if (.not. avg) text_lname=reminpocsm_long_name
        if (      avg) text_lname='averaged '// reminpocsm_long_name
        varid =
     &  nccreate(ncid,reminpocsm_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',reminpocsm_units)
      endif


      !remin of DOC
      if ( ( wrt_remindoc .and. .not. avg) .or. (wrt_avg_remindoc .and. avg) ) then
        if (.not. avg) text_lname=remindoc_long_name
        if (      avg) text_lname='averaged '// remindoc_long_name
        varid =
     &  nccreate(ncid,remindoc_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',remindoc_units)
      endif


      !remin of DOC refractory
      if ( ( wrt_remindocref .and. .not. avg) .or. (wrt_avg_remindocref .and. avg) ) then
        if (.not. avg) text_lname=remindocref_long_name
        if (      avg) text_lname='averaged '// remindocref_long_name
        varid =
     &  nccreate(ncid,remindocref_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',remindocref_units)
      endif


      end subroutine def_vars_mag  !]
      ! ----------------------------------------------------------------------

     ! ----------------------------------------------------------------------
      subroutine calc_avg_mag ![
      implicit none

      ! local
      real    :: coef

      if (wrt_file_avg) then 

      navg_mag = navg_mag +1
      coef = 1./navg_mag

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'mag :: started averaging. ',
     &   'output_period_avg (s) =', output_period_avg
      endif

      t_avg_mag = t_avg_mag*(1-coef) + time*coef


      !2d biomass
      if (wrt_avg_b2d) then 
         algae_b2d_avg(i0:i1,j0:j1) = algae_b2d_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + algae_b2d(i0:i1,j0:j1) * coef
      endif

      !3d biomass
      if (wrt_avg_b) then 
         algae_b_avg(i0:i1,j0:j1,:) = algae_b_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_b(i0:i1,j0:j1,:) * coef
      endif

      !leaf-area density
      if (wrt_avg_lad) then 
         algae_lad_avg(i0:i1,j0:j1,:) = algae_lad_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_lad(i0:i1,j0:j1,:) * coef
      endif

 
      !length
      if (wrt_avg_length) then 
         algae_length_avg(i0:i1,j0:j1) = algae_length_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + algae_length(i0:i1,j0:j1) * coef
      endif

      !shape function
      if (wrt_avg_shape) then 
         algae_shape_avg(i0:i1,j0:j1,:) = algae_shape_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_shape(i0:i1,j0:j1,:) * coef
      endif

      !npp
      if (wrt_avg_npp) then 
         algae_npp_avg(i0:i1,j0:j1) = algae_npp_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + algae_npp(i0:i1,j0:j1) * coef
      endif

      !loss
      if (wrt_avg_loss) then 
         algae_loss_avg(i0:i1,j0:j1) = algae_loss_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + algae_loss(i0:i1,j0:j1) * coef
      endif


      !harvest
      if (wrt_avg_harv) then 
         algae_harv_avg(i0:i1,j0:j1) = algae_harv_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + algae_harv(i0:i1,j0:j1) * coef
      endif

      !nitrate uptake
      if (wrt_avg_uptno3) then 
         algae_uptno3_avg(i0:i1,j0:j1,:) = algae_uptno3_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_uptno3(i0:i1,j0:j1,:) * coef
      endif


      !exudation loss
      if (wrt_avg_exuL) then 
         algae_exuL_avg(i0:i1,j0:j1,:) = algae_exuL_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_exuL(i0:i1,j0:j1,:) * coef
      endif

      !small particulate loss
      if (wrt_avg_spomL) then 
         algae_spomL_avg(i0:i1,j0:j1,:) = algae_spomL_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_spomL(i0:i1,j0:j1,:) * coef
      endif


      !frond particulate loss
      if (wrt_avg_frondL) then 
         algae_frondL_avg(i0:i1,j0:j1,:) = algae_frondL_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_frondL(i0:i1,j0:j1,:) * coef
      endif


      !3D PAR
      if (wrt_avg_par) then 
         PAR_3d_mag_avg(i0:i1,j0:j1,:) = PAR_3d_mag_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + PAR_3d_mag(i0:i1,j0:j1,:) * coef
      endif

      ! remin of frond poc 
      if (wrt_avg_reminpoclrg) then 
         remin_lrgpocalgae_avg(i0:i1,j0:j1,:) = remin_lrgpocalgae_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + remin_lrgpocalgae(i0:i1,j0:j1,:) * coef
      endif

      ! remin of small poc 
      if (wrt_avg_reminpocsm) then 
         remin_smpocalgae_avg(i0:i1,j0:j1,:) = remin_smpocalgae_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + remin_smpocalgae(i0:i1,j0:j1,:) * coef
      endif


      ! remin of doc 
      if (wrt_avg_remindoc) then 
         remin_docalgae_avg(i0:i1,j0:j1,:) = remin_docalgae_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + remin_docalgae(i0:i1,j0:j1,:) * coef
      endif

      ! remin of doc refractory 
      if (wrt_avg_remindocref) then 
         remin_docrefalgae_avg(i0:i1,j0:j1,:) = remin_docrefalgae_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + remin_docrefalgae(i0:i1,j0:j1,:) * coef
      endif



      endif !wrt_file_avg

      end subroutine calc_avg_mag !]
      !-----------------------------------------------------------------

      !------------- OUT DATED SUBROUTINE ----------------
      subroutine conv_mag_frc
      ! =========================================
      !  Convert mag forcing (biomass, sflx_pockelp)
      ! to what is needed for ROMS-BEC coupling 


      !                What this subroutine does
      !    1. converts 2D biomass to 3D biomass (computes algae_b)
      !          - computes algae_shape (shape function)
      !          - algae_b = algae_b2d * algae_shape
      !    2. computes leaf-area density
      !          - computes algae_lad, which needs algae_length
      !    3. converts surface flux of sinking biomass to BEC-compatible
      !       units
      ! =========================================
      implicit none
      ! inputs
      integer :: i,j,k
      real :: conv


      ! presumably can just work on interior points, since no algae
      ! variables will be placed on stencil that requires buffers
      ! (...but maybe for drag???!?)
       do j=1,ny
         do i=1,nx
            !Only compute kelp characteristics if there is biomass there
            if (algae_b2d(i,j).ne.0) then     
               
               !----------------------------------------------------
               !Define cultivation sigma level based on F_z_cult
               !(defined by user in mag.opt)
               !----------------------------------------------------
               !call find_k(i,j,algae_zbot*-1.0,kbot(i,j))
               !call find_k(i,j,-algae_zbot,kbot)
               !print *, 'kbot =', kbot 
               !        Print to chek that find_k works--> looks good
               !        (9/13/2024 9AM)
               !print *, 'algae_zbot=', algae_zbot              
               !print *,'kbot=',  kbot
               !print *, 'zr(kbot)=', z_r(i,j,kbot)

               !----------------------
               !Calculate algae_length (based on z-integrated biomass)
               !-----------------------
               call calc_algaelength(i,j)

               ! -----------------------------------------------------
               ! Convert vert. integral of biomass to a vert. profile
               ! ----------------------------------------------------
               ! computes algae_shape(i,j)
               call calc_shape(i,j)
               !Now distribute biomass
               do k=1,N
                  algae_b(i,j,k) = algae_b2d(i,j) * algae_shape(i,j,k)
               enddo
               

               ! ---------------------------------------
               ! Convert biomass(z) to leaf-area density
               ! ---------------------------------------
               call calc_LAD(i,j)

               !---------------------------------
               ! Convert biomass to fixed nitrogen
               ! Nf used for PAR modulation
               !----------------------------------
               !print *, 'algae_Nf=', algae_Nf(i,j,kbot+1)
               algae_Nf(i,j,:) = Qmin * algae_b(i,j,:)

               
               !----------------------------------------
               !Compute exudation rate in nitrogen units
               !----------------------------------------
               !    mmol N / m^3 / sec      mg N / m^3      g/mol   1/s
               !algae_em_nitrogen(i,j,:) = algae_Nf(i,j,:) * MW_N * exu0


            endif

            ! --------------------------------------
            !         Convert pockelp surface flux to 
            !         units of [mmol C / m^3 / sec]
            !         from     [g C / m^2 / day]
            ! --------------------------------------
            !Convert g C / m^2 / day to mmol C / m^2 / sec
                   ! g C --> mmol C   ! 1/day --> 1/sec
            conv = (1000. / 12.01) * (1./86400.)
            !divide by Hz at surface to make it a source term at surface and
            !give units of per m^3
            !print *, 'Hz(i,j,N) = ', Hz(i,j,N)
            pockelp_prod_surf(i,j) = sflx_pockelp_mass(i,j)*conv / Hz(i,j,N)

        enddo !< -- i
      enddo !<-- j

      end subroutine conv_mag_frc 



#endif /* MAG for whole module */
      end module mag  
