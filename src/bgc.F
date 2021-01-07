      module bgc

      ! BGC - bio-geo-chemical module
      ! -----------------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Dec)
      ! Refactoring of ETH Zurich roms code with BEC2, code which has been used
      ! by Pierre Damien to run full pacific models.

#include "cppdefs.h"
#ifdef BIOLOGY_BEC2

      use read_write
      use tracers ! For iPO4, etc, indices of bgc tracers (formerly in param.h)

      implicit none

      private ! Make all variable private to this module unless public specified

      ! Includes:
      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"
      ! ncvars.h for nchis, ncavg, ncrst
#include "ncvars.h"
#include "bgc_ecosys_bec2.h"

      ! ****************************************************************
      ! User inputs

      integer, parameter :: nbgc_flx = 2 ! Number of surface bgc flux forcings
#ifdef PCO2AIR_FORCING
     &                                +1
      ! pco2air forcing
      character(len=8) :: pco2air_vname = 'pco2_air'
      character(len=9) :: pco2air_tname = 'pco2_time'
#endif

      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_frc = 0 ! Factor 2 only for now
      ! dust forcing
      character(len=4) :: dust_vname = 'dust'
      character(len=9) :: dust_tname = 'dust_time'
      ! iron forcing
      character(len=4) :: iron_vname = 'iron'
      character(len=9) :: iron_tname = 'iron_time'

#include "bgc_tracers.h"
#include "bgc_bio_diag.h"

      ! End of user inputs
      ! *************************************************************

      ! Netcdf variables for surface forcing - BGC FLUX ONLY
      real                         :: flx_data(GLOBAL_2D_ARRAY,2,nbgc_flx)
      integer, dimension(nbgc_flx) :: flx_grd_type  = 0
      integer, dimension(nbgc_flx) :: flx_file_indx = 0
      integer, dimension(nbgc_flx) :: flx_irec      = 0
      integer, dimension(nbgc_flx) :: flx_it1 = 1, flx_it2 = 2
      real, dimension(2,nbgc_flx)  :: flx_times = -99 ! [-99,-99]

      ! Misc from porting ETH code:
      character(len=40) tclm_name(NT) ! Taken from ETH code ncvars.h

      public set_bgc_surf_frc
#ifdef PCO2AIR_FORCING
      public set_pco2air_tile
#endif
#ifdef BEC2_DIAG
      public def_his_bgc_diag
      public def_rst_bgc_diag
      public wrt_his_bgc_diag
      public wrt_rst_bgc_diag
# ifdef AVERAGES
      public def_avg_bgc_diag
      public wrt_avg_bgc_diag
      public set_avg_bgc_diag
# endif
#endif
      public init_scalars_bec2

      contains
! ----------------------------------------------------------------------

      subroutine set_bgc_surf_frc(istr,iend,jstr,jend)
      ! read in bgc surface flux and interpolate to model time.
      ! Taken from get_smth.F & set_smth.F of ETH code.

      ! Since bgc is set up such that iron and dust have flux at surface but
      ! are not directly bgc tracers, I have separated them into set_bgc_frc

      use read_write, only: set_frc_var_tile

      implicit none

! Needed for iron and dust variables, not yet ported to this module. Should do.
#include "bgc_forces.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! dust & iron controlled by # if defined BIOLOGY_BEC2 in set_forces.F of ETH code

      ! dust
      call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        dust_vname,         dust_tname,       ! Text name of array and time
     &        flx_data(:,:,:,1),  dust,             ! Variable arrays
     &        flx_times(:,1),     flx_grd_type(1),  ! Input times & Grid type
     &        flx_file_indx(1),   flx_irec(1),      ! File indx & Current input record
     &        flx_it1(1),         flx_it2(1),       ! Time index placeholders
     &        interp_frc )                          ! Online spatial interpolation flag

      ! iron
      call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        iron_vname,         iron_tname,       ! Text name of array and time
     &        flx_data(:,:,:,2),  iron,             ! Variable arrays
     &        flx_times(:,2),     flx_grd_type(2),  ! Input times & Grid type
     &        flx_file_indx(2),   flx_irec(2),      ! File indx & Current input record
     &        flx_it1(2),         flx_it2(2),       ! Time index placeholders
     &        interp_frc )                          ! Online spatial interpolation flag

!      if(mynode==0) print*, 'dust(20,20) ',dust(20,20)
!      if(mynode==0) print*, 'iron(20,20) ',iron(20,20)

      end subroutine set_bgc_surf_frc
! ----------------------------------------------------------------------

#ifdef PCO2AIR_FORCING
      subroutine set_pco2air_tile(istr,iend,jstr,jend)
      ! read in bgc pco2air flux and interpolate to model time.
      ! Taken from get_forces.F, get_smth.F & set_smth.F of ETH code.

      use read_write, only: set_frc_var_tile

      implicit none

! Needed for pco2air variable, not yet ported to this module. Should do so at some point.
#include "bgc_forces.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! pco2air controlled by # if defined BIOLOGY_BEC2 in set_forces.F of ETH code

      ! Index 3 below is because currently there are 3 surface forcings,
      ! and dust=1 and iron=2.
      call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        pco2air_vname,      pco2air_tname,    ! Text name of array and time
     &        flx_data(:,:,:,3),  pco2air,          ! Variable arrays
     &        flx_times(:,3),     flx_grd_type(3),  ! Input times & Grid type
     &        flx_file_indx(3),   flx_irec(3),      ! File indx & Current input record
     &        flx_it1(3),         flx_it2(3),       ! Time index placeholders
     &        interp_frc )                          ! Online spatial interpolation flag

!      if(mynode==0) print*, 'pco2air(20,20) ',pco2air(20,20)

      end subroutine set_pco2air_tile
#endif /* PCO2AIR_FORCING */
! ----------------------------------------------------------------------

#ifdef BEC2_DIAG
      subroutine def_his_bgc_diag(r2dgrd
#ifdef SOLVE3D
     &                           ,r3dgrd
#endif
     &                                  )

      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F

      implicit none

      ! Inputs
      integer r2dgrd(3) !, u2dgrd(3), v2dgrd(3)
#ifdef SOLVE3D
     &      , r3dgrd(4) !, u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif
      ! Local
      integer idiag, ierr

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrtavg_bec2_diag_2d(idiag)) then
          call nc_define_var(nchis, vname_bec2_diag_2d(1,idiag),
     &    vname_bec2_diag_2d(2,idiag), vname_bec2_diag_2d(3,idiag),
     &    r2dgrd, ierr )
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrtavg_bec2_diag_3d(idiag)) then
          call nc_define_var(nchis, vname_bec2_diag_3d(1,idiag),
     &    vname_bec2_diag_3d(2,idiag), vname_bec2_diag_3d(3,idiag),
     &    r3dgrd, ierr )
        endif
      enddo

      end subroutine def_his_bgc_diag
#endif /* BEC2_DIAG */
! ----------------------------------------------------------------------

! Need AVERAGES as 'ncavg' only exists in AVERAGES is defined.
#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine def_avg_bgc_diag(r2dgrd
# ifdef SOLVE3D
     &                           ,r3dgrd
# endif
     &                                  )

      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F

      implicit none

      ! Inputs
      integer r2dgrd(3) !, u2dgrd(3), v2dgrd(3)
# ifdef SOLVE3D
     &      , r3dgrd(4) !, u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif
      ! Local
      integer idiag, ierr
      character(len=64) long_name

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrtavg_bec2_diag_2d(idiag)) then
          long_name='averaged '/ /vname_bec2_diag_2d(2,idiag) ! Add averaged to long name
          call nc_define_var(ncavg, vname_bec2_diag_2d(1,idiag),
     &    long_name, vname_bec2_diag_2d(3,idiag), r2dgrd, ierr )
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrtavg_bec2_diag_3d(idiag)) then
          long_name='averaged '/ /vname_bec2_diag_3d(2,idiag) ! Add averaged to long name
          call nc_define_var(ncavg, vname_bec2_diag_3d(1,idiag),
     &    long_name, vname_bec2_diag_3d(3,idiag), r3dgrd, ierr )
        endif
      enddo

      end subroutine def_avg_bgc_diag
#endif /* BEC2_DIAG */
#endif /* AVERAGES */
! ----------------------------------------------------------------------

#ifdef BEC2_DIAG
      subroutine def_rst_bgc_diag(r2dgrd
#ifdef SOLVE3D
     &                           ,r3dgrd
#endif
     &                                  )

      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F

      implicit none

      ! Inputs
      integer r2dgrd(3) !, u2dgrd(3), v2dgrd(3)
#ifdef SOLVE3D
     &      , r3dgrd(4) !, u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif
      ! Local
      integer idiag, ierr

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrtavg_bec2_diag_2d(idiag)) then
          call nc_define_var(ncrst, vname_bec2_diag_2d(1,idiag),
     &      vname_bec2_diag_2d(2,idiag), vname_bec2_diag_2d(3,idiag),
     &      r2dgrd, ierr )
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrtavg_bec2_diag_3d(idiag)) then
          call nc_define_var(ncrst, vname_bec2_diag_3d(1,idiag),
     &      vname_bec2_diag_3d(2,idiag), vname_bec2_diag_3d(3,idiag),
     &      r3dgrd, ierr )
        endif
      enddo

      end subroutine def_rst_bgc_diag
#endif /* BEC2_DIAG */
! ----------------------------------------------------------------------

#ifdef BEC2_DIAG
      subroutine wrt_his_bgc_diag (record, ierr)
      ! Write BEC2_DIAG to history file
      implicit none

      ! inputs
      integer, intent(in) :: record ! The timestep to be recorded
      integer :: ierr ! Error handling
      ! local
      integer itrc

      ! 2d diagnostics:
      do itrc=1,nr_bec2_diag_2d
         if (wrtavg_bec2_diag_2d(itrc)) then
           call nc_write_var(nchis, bec2_diag_2d(START_2D_ARRAY,itrc), 1,
     &             vname_bec2_diag_2d(1,itrc), rp_var, record, ierr )
         endif
      end do

      ! 3d diagnostics:
      do itrc=1,nr_bec2_diag_3d
         if (wrtavg_bec2_diag_3d(itrc)) then
           call nc_write_var(nchis, bec2_diag_3d(START_2D_ARRAY,1,itrc), N,
     &             vname_bec2_diag_3d(1,itrc), rp_var, record, ierr )
         endif
      end do

      end subroutine wrt_his_bgc_diag
#endif /* BEC2_DIAG */
! ----------------------------------------------------------------------

#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine wrt_avg_bgc_diag (record, ierr)
      ! Write BEC2_DIAG to avg file
      implicit none

      ! inputs
      integer, intent(in) :: record ! The timestep to be recorded
      integer :: ierr ! Error handling
      ! local
      integer itrc

      ! 2d diagnostics:
      do itrc=1,nr_bec2_diag_2d
         if (wrtavg_bec2_diag_2d(itrc)) then
           call nc_write_var(ncavg, bec2_diag_2d_avg(START_2D_ARRAY,itrc),
     &          1, vname_bec2_diag_2d(1,itrc), rp_var, record, ierr )
         endif
      end do

      ! 3d diagnostics:
      do itrc=1,nr_bec2_diag_3d
         if (wrtavg_bec2_diag_3d(itrc)) then
           call nc_write_var(ncavg, bec2_diag_3d_avg(START_2D_ARRAY,1,itrc),
     &          N, vname_bec2_diag_3d(1,itrc), rp_var, record, ierr )
         endif
      end do

      end subroutine wrt_avg_bgc_diag
#endif /* BEC2_DIAG */
#endif /* AVERAGES */
! ----------------------------------------------------------------------

#ifdef BEC2_DIAG
      subroutine wrt_rst_bgc_diag (record, ierr)
      ! Write BEC2_DIAG to restart file
      implicit none

      ! inputs
      integer, intent(in) :: record ! The timestep to be recorded
      integer :: ierr ! Error handling
      ! local
      integer itrc

      ! 2d diagnostics:
      do itrc=1,nr_bec2_diag_2d
        call nc_write_var(ncrst, bec2_diag_2d(START_2D_ARRAY,itrc),
     &          1, vname_bec2_diag_2d(1,itrc), rp_var, record, ierr )
      end do

      ! 3d diagnostics:
      do itrc=1,nr_bec2_diag_3d
        call nc_write_var(ncrst, bec2_diag_3d(START_2D_ARRAY,1,itrc),
     &          N, vname_bec2_diag_3d(1,itrc), rp_var, record, ierr )
      end do

      end subroutine wrt_rst_bgc_diag
#endif /* BEC2_DIAG */
! ----------------------------------------------------------------------

#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine set_avg_bgc_diag(istrR,iendR,jstrR,jendR,cff,cff1)

      implicit none

      ! inputs
      integer istrR,iendR,jstrR,jendR
      real cff,cff1
      ! local
      integer i,j,k,itrc

      do itrc=1,nr_bec2_diag_2d
         do j=jstrR,jendR
            do i=istrR,iendR
               bec2_diag_2d_avg(i,j,itrc) = cff1*bec2_diag_2d_avg(i,j,itrc)
     &                                    + cff*bec2_diag_2d(i,j,itrc)
            end do
         end do
      end do
      do itrc=1,nr_bec2_diag_3d
        do k=1,N
          do j=jstrR,jendR
             do i=istrR,iendR
                bec2_diag_3d_avg(i,j,k,itrc) = cff1*bec2_diag_3d_avg(i,j,k,itrc)
     &                                       + cff*bec2_diag_3d(i,j,k,itrc)
             end do
          end do
        end do
      end do

      end subroutine set_avg_bgc_diag
#endif /* BEC2_DIAG */
#endif /* AVERAGES */
! ----------------------------------------------------------------------

      subroutine init_scalars_bec2()
!
! Set initial values for  globally accessible (stored in common
! blocks) scalar variables of the BEC model.
!
      implicit none
      integer i,j, itrc ,lvar,lenstr ! ierr, DevinD removed ierr as arguement
      integer omp_get_num_threads
!#include "param.h"
#include "scalars.h"
!#include "ncvars.h"
#include "comm_vars.h"
#include "ncroms_bgcflux.h"

! CN: take this part out of the NEVER_DEF part
!#ifdef BEC2_DIAG
!# include "ecosys_bec2.h"
!#endif

! DevinD commented out 20201215
!#ifdef NEVER_DEF
!
!!
!! Find out how many threads are created and check whether the number
!! of threads exceeds maximum allowed, as well as check if the number
!! of tiles is divisible by the number of threads, so the job can be
!! evenly partitioned, complain and signal to terminate, if something
!! is wrong.
!!
!C$OMP PARALLEL
!C$OMP CRITICAL (isca_cr_rgn)
!      numthreads=1
!C$      numthreads=omp_get_num_threads()
!C$OMP END CRITICAL (isca_cr_rgn)
!C$OMP END PARALLEL
!      MPI_master_only write(*,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
!     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
!
!#ifdef MPI
!# define NUMPROC NNODES
!#else
!# define NUMPROC 1
!#endif
!
!      if (numthreads.gt.NUMPROC) then
!        MPI_master_only write(*,'(/1x,A,I3/)')
!     &    'ERROR: Requested number of threads exceeds setting: ',NUMPROC
!        ierr=ierr+1                                      !--> ERROR
!      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
!        MPI_master_only write(*,
!     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
!     &                'wrong choice of numthreads =', numthreads,
!     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
!        ierr=ierr+1                                      !--> ERROR
!      endif
!!
!      do i=0,2*NT+1
!        global_sum(i)=QuadZero
!      enddo
!#endif /* NEVER_DEF */

!
! Define names of variables in NetCDF output files.
! Inner dimension is for variable type and outer is for variable
! attributes (name, long-name, units, field).
!

      ! DevinD: NOTE: 'vname' was previously used to store all
      ! variable names (primary vars, tracers, etc).
      ! However, tracers variable names, etc are now stored in tracers.F
      ! in t_vname, t_lname, etc. Hence populating vname for bgc should not
      ! be necessary. I have left it populated for now as it is somewhat harmless.
      ! HOWEVER, vname_bec2_diag_2d & vname_bec2_diag_3d are STILL used and must
      ! be correctly populated.

#ifdef SOLVE3D
      vname(1,indxPO4)='PO4                                       '
      vname(2,indxpO4)='Phosphate                                 '
      vname(3,indxpO4)='mMol P m-3                                '
      vname(4,indxpO4)='  '

CR      write(*,'(1x,A,I3)') 'indxNO3 =',indxNO3
      vname(1,indxNO3)='NO3                                       '
      vname(2,indxNO3)='Nitrate                                   '
      vname(3,indxNO3)='mMol N m-3                                '
      vname(4,indxNO3)='  '

      vname(1,indxSIO3)='SiO3                                      '
      vname(2,indxSIO3)='Silicate                                  '
      vname(3,indxSIO3)='mMol Si m-3                               '
      vname(4,indxSIO3)='  '

CR      write(*,'(1x,A,I3)') 'indxNH4 =',indxNH4
      vname(1,indxNH4)='NH4                                       '
      vname(2,indxNH4)='Ammonium                                  '
      vname(3,indxNH4)='mMol N m-3                                '
      vname(4,indxNH4)='  '


      vname(1,indxFE)='Fe                                        '
      vname(2,indxFE)='Iron                                      '
      vname(3,indxFE)='mMol Fe m-3                               '
      vname(4,indxFE)='  '

      vname(1,indxO2)='O2                                        '
      vname(2,indxO2)='Oxygen                                    '
      vname(3,indxO2)='mMol O2 m-3                               '
      vname(4,indxO2)='  '

      vname(1,indxDIC)='DIC                                       '
      vname(2,indxDIC)='Dissolved inorganic carbon                '
      vname(3,indxDIC)='mMol C m-3                                '
      vname(4,indxDIC)='  '

      vname(1,indxALK)='Alk                                       '
      vname(2,indxALK)='Alkalinity                                '
      vname(3,indxALK)='mMol m-3                                '
      vname(4,indxALK)='  '

      vname(1,indxDOC)='DOC                                       '
      vname(2,indxDOC)='Dissolved organic carbon                  '
      vname(3,indxDOC)='mMol C m-3                                '
      vname(4,indxDOC)='  '

      vname(1,indxDon)='DON                                      '
      vname(2,indxDon)='Dissolved organic nitrogen               '
      vname(3,indxDon)='mMol N m-3                               '
      vname(4,indxDon)='  '

      vname(1,indxDofe)='DOFE                                     '
      vname(2,indxDofe)='Dissolved organic iron                   '
      vname(3,indxDofe)='mMol Fe m-3                              '
      vname(4,indxDofe)='  '

      vname(1,indxDop)='DOP                                      '
      vname(2,indxDop)='Dissolved organic phosphorus             '
      vname(3,indxDop)='mMol P m-3                               '
      vname(4,indxDop)='  '

      vname(1,indxDopr)='DOPR'
      vname(2,indxDopr)='Refractory dissolved organic phosphorus  '
      vname(3,indxDopr)='mMol P m-3'
      vname(4,indxDopr)='  '

      vname(1,indxDonr)='DONR'
      vname(2,indxDonr)='Refractory dissolved organic nitrogen    '
      vname(3,indxDonr)='mMol N m-3'
      vname(4,indxDonr)='  '

      vname(1,indxZOOC)='ZOOC'
      vname(2,indxZOOC)='Zooplankton'
      vname(3,indxZOOC)='mMol C m-3'
      vname(4,indxZOOC)='  '

      vname(1,indxSPCHL)='SPCHL'
      vname(2,indxSPCHL)='Small phytoplankton chlorophyll'
      vname(3,indxSPCHL)='mg Chl-a m-3'
      vname(4,indxSPCHL)='  '

      vname(1,indxSPC)='SPC                                       '
      vname(2,indxSPC)='Small phytoplankton carbon                '
      vname(3,indxSPC)='mMol C m-3                                '
      vname(4,indxSPC)='  '

      vname(1,indxSPFE)='SPFE                                     '
      vname(2,indxSPFE)='Small phytoplankton iron                 '
      vname(3,indxSPFE)='mMol Fe m-3                              '
      vname(4,indxSPFE)='  '

      vname(1,indxSPCACO3)='SPCACO3                                   '
      vname(2,indxSPCACO3)='Small phytoplankton CaCO3                 '
      vname(3,indxSPCACO3)='mMol CaCO3 m-3                            '
      vname(4,indxSPCACO3)='  '

      vname(1,indxDIATCHL)='DIATCHL                                   '
      vname(2,indxDIATCHL)='Diatom chlorophyll                        '
      vname(3,indxDIATCHL)='mg Chl-a m-3                              '
      vname(4,indxDIATCHL)='  '

      vname(1,indxDIATC)='DIATC                                     '
      vname(2,indxDIATC)='Diatom carbon                             '
      vname(3,indxDIATC)='mMol C m-3                                '
      vname(4,indxDIATC)='  '

      vname(1,indxDIATFE)='DIATFE                                    '
      vname(2,indxDIATFE)='Diatom Iron                               '
      vname(3,indxDIATFE)='mMol Fe m-3                               '
      vname(4,indxDIATFE)='  '

      vname(1,indxDIATSI)='DIATSI                                    '
      vname(2,indxDIATSI)='Diatom silicon                            '
      vname(3,indxDIATSI)='mMol Si m-3                               '
      vname(4,indxDIATSI)='  '

      vname(1,indxDiazchl)='DIAZCHL                                  '
      vname(2,indxDiazchl)='Diazotroph chlorophyll                   '
      vname(3,indxDiazchl)='mg Chl-a m-3                             '
      vname(4,indxDiazchl)='  '

      vname(1,indxdiazc)='DIAZC                                    '
      vname(2,indxDiazc)='Diazotroph carbon                        '
      vname(3,indxDiazc)='mMol C m-3                               '
      vname(4,indxDiazc)='  '

      vname(1,indxDiazfe)='DIAZFE                                   '
      vname(2,indxDiazfe)='Diazotroph iron                          '
      vname(3,indxDiazfe)='mMol Fe m-3                              '
      vname(4,indxDiazfe)='  '

# ifdef BEC_COCCO
      vname(1,indxCoccoc)='COCCOC                                   '
      vname(2,indxCoccoc)='Coccolithophore Carbon                   '
      vname(3,indxCoccoc)='mMol C m-3                              '
      vname(4,indxCoccoc)='  '

      vname(1,indxCoccochl)='COCCOCHL                                 '
      vname(2,indxCoccochl)='Coccolithophore Chlorophyll              '
      vname(3,indxCoccochl)='mg Chl-a m-3                              '
      vname(4,indxCoccochl)='  '

      vname(1,indxCoccocal)='COCCOCAL                                '
      vname(2,indxCoccocal)='Coccolithophore CaCO3                   '
      vname(3,indxCoccocal)='mMol CaCO3 m-3                          '
      vname(4,indxCoccocal)='  '

      vname(1,indxCoccofe)='COCCOFE                                   '
      vname(2,indxCoccofe)='Coccolithophore iron                    '
      vname(3,indxCoccofe)='mMol Fe m-3                              '
      vname(4,indxCoccofe)='  '

      vname(1,indxCal)='CAL                                   '
      vname(2,indxCal)='CaCO3                                   '
      vname(3,indxCal)='mMol CaCO3 m-3                         '
      vname(4,indxCal)='  '
# endif

# ifdef Ncycle_SY
      vname(1,indxno2)='NO2                                   '
      vname(2,indxno2)='Nitrite                                 '
      vname(3,indxno2)='mMol N m-3                              '
      vname(4,indxno2)='  '

      vname(1,indxn2)='N2                              '
      vname(2,indxn2)='Dinitrogen   '
      vname(3,indxn2)='mMol N2 m-3                         '
      vname(4,indxn2)='  '

      vname(1,indxn2o)='N2O                              '
      vname(2,indxn2o)='Nitrous oxide   '
      vname(3,indxn2o)='mMol N2O m-3                         '
      vname(4,indxn2o)='  '

# ifdef N2O_TRACER_DECOMP
      vname(1,indxn2o_ao1)='N2O_AO1                                 '
      vname(2,indxn2o_ao1)='Nitrous oxide produced via AO1             '
      vname(3,indxn2o_ao1)='mMol N2O m-3                              '
      vname(4,indxn2o_ao1)='  '

      vname(1,indxn2o_siden)='N2O_SIDEN                                '
      vname(2,indxn2o_siden)='Nitrous oxide consumed via denitrif. '
      vname(3,indxn2o_siden)='mMol N2O m-3                             '
      vname(4,indxn2o_siden)='  '

      vname(1,indxn2o_soden)='N2O_SODEN                              '
      vname(2,indxn2o_soden)='Nitrous oxide produced via denitrif.   '
      vname(3,indxn2o_soden)='mMol N2O m-3                         '
      vname(4,indxn2o_soden)='  '

      vname(1,indxn2o_atm)='N2O_ATM                              '
      vname(2,indxn2o_atm)='N2O of atmospheric provenance   '
      vname(3,indxn2o_atm)='mMol N2 m-3                         '
      vname(4,indxn2o_atm)='  '

      vname(1,indxn2_sed)='N2_SED                              '
      vname(2,indxn2_sed)='N2 from sedimentary denitrification   '
      vname(3,indxn2_sed)='mMol N2 m-3                         '
      vname(4,indxn2_sed)='  '
# endif
# endif

# ifdef N2O_NEV
      vname(1,indxn2o_nev)='N2O_NEV                              '
      vname(2,indxn2o_nev)='N2O (Nevison)   '
      vname(3,indxn2o_nev)='mMol N2O m-3                         '
      vname(4,indxn2o_nev)='  '
# endif

      ! DevinD: commented these out as were overidding Akv & Akt vname's
!      vname(1,indxdust)='dust                                       '
!      vname(2,indxdust)='dust flux                      '
!      vname(3,indxdust)='kg m-2 s-1                                '
!      vname(4,indxdust)='  '
!
!      vname(1,indxiron)='iron                                    '
!      vname(2,indxiron)='iron flux                      '
!      vname(3,indxiron)='nmol cm-2 s-1                             '
!      vname(4,indxiron)='  '

# ifndef BEC2_DIAG
      vname(1,indxPH)='PH'
      vname(2,indxPH)='surf pH value'
      vname(3,indxPH)=' '
      vname(4,indxPH)=' '
      vname(1,indxPCO2)='pCO2                                 '
      vname(2,indxPCO2)='Surface water pCO2                   '
      vname(3,indxPCO2)='ppm                                  '
      vname(4,indxPCO2)='  '
      vname(1,indxPCO2air)='pCO2air                           '
      vname(2,indxPCO2air)='Atmospheric pCO2                  '
      vname(3,indxPCO2air)='ppm                               '
      vname(4,indxPCO2air)='           '
      vname(1,indxPARinc)='PARinc                                '
      vname(2,indxPARinc)=
     &  'Incoming photosynth. active radiation '
      vname(3,indxPARinc)='W m-2                                 '
      vname(4,indxPARinc)='  '
      vname(1,indxPAR)='PAR'
      vname(2,indxPAR)=
     &  'Photosynthetically active radiation'
      vname(3,indxPAR)='W m-2'
      vname(4,indxPAR)='  '
# endif /*BEC2_DIAG */



# ifdef BEC2_DIAG
      !
      ! 2d diagnotics:
      !
#  if defined CCHEM_MOCSY && defined CCHEM_TODEPTH
      vname_bec2_diag_3d(1,indxPH)='PH'
      vname_bec2_diag_3d(2,indxPH)='surf pH value'
      vname_bec2_diag_3d(3,indxPH)=' '
      vname_bec2_diag_3d(4,indxPH)=' '
      vname_bec2_diag_3d(1,indxPCO2)='pCO2'
      vname_bec2_diag_3d(2,indxPCO2)='Surface water pCO2'
      vname_bec2_diag_3d(3,indxPCO2)='ppm'
      vname_bec2_diag_3d(4,indxPCO2)='  '
#  else
      vname_bec2_diag_2d(1,indxPH)='PH'
      vname_bec2_diag_2d(2,indxPH)='surf pH value'
      vname_bec2_diag_2d(3,indxPH)=' '
      vname_bec2_diag_2d(4,indxPH)=' '
      vname_bec2_diag_2d(1,indxPCO2)='pCO2'
      vname_bec2_diag_2d(2,indxPCO2)='Surface water pCO2'
      vname_bec2_diag_2d(3,indxPCO2)='ppm'
      vname_bec2_diag_2d(4,indxPCO2)='  '
#  endif
      vname_bec2_diag_2d(1,indxPCO2air)='pCO2air'
      vname_bec2_diag_2d(2,indxPCO2air)='Atmospheric pCO2'
      vname_bec2_diag_2d(3,indxPCO2air)='ppm'
      vname_bec2_diag_2d(4,indxPCO2air)='           '
      vname_bec2_diag_2d(1,indxPARinc)='PARinc'
      vname_bec2_diag_2d(2,indxPARinc)=
     &  'Incoming photosynth. active radiation '
      vname_bec2_diag_2d(3,indxPARinc)='W m-2'
      vname_bec2_diag_2d(4,indxPARinc)='  '
      vname_bec2_diag_2d(1,indxFGO2)='FG_O2'
      vname_bec2_diag_2d(2,indxFGO2)='Air-sea flux of O2'
      vname_bec2_diag_2d(3,indxFGO2)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFGO2)=' '
      vname_bec2_diag_2d(1,indxFGCO2)='FG_CO2'
      vname_bec2_diag_2d(2,indxFGCO2)='Air-sea flux of CO2'
      vname_bec2_diag_2d(3,indxFGCO2)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxFGCO2)=' '
      vname_bec2_diag_2d(1,indxWS10m)='WS10m'
      vname_bec2_diag_2d(2,indxWS10m)='Wind speed at 10m'
      vname_bec2_diag_2d(3,indxWS10m)='m/s '
      vname_bec2_diag_2d(4,indxWS10m)=' '
      vname_bec2_diag_2d(1,indxXKW)='XKW'
      vname_bec2_diag_2d(2,indxXKW)='XKW'
      vname_bec2_diag_2d(3,indxXKW)='m/s '
      vname_bec2_diag_2d(4,indxXKW)=' '
      vname_bec2_diag_2d(1,indxATMPRESS)='ATM_PRESS'
      vname_bec2_diag_2d(2,indxATMPRESS)='Atmospheric pressure'
      vname_bec2_diag_2d(3,indxATMPRESS)='atm '
      vname_bec2_diag_2d(4,indxATMPRESS)=' '

      vname_bec2_diag_2d(1,indxSCHMIDTO2)='SCHMIDT_O2'
      vname_bec2_diag_2d(2,indxSCHMIDTO2)='Schmidt number for O2'
      vname_bec2_diag_2d(3,indxSCHMIDTO2)=' '
      vname_bec2_diag_2d(4,indxSCHMIDTO2)=' '

      vname_bec2_diag_2d(1,indxO2SAT)='O2SAT'
      vname_bec2_diag_2d(2,indxO2SAT)='O2 saturation concentration'
      vname_bec2_diag_2d(3,indxO2SAT)='mmol/m3 '
      vname_bec2_diag_2d(4,indxO2SAT)=' '
      vname_bec2_diag_2d(1,indxSCHMIDTCO2)='SCHMIDT_CO2'
      vname_bec2_diag_2d(2,indxSCHMIDTCO2)='Schmidt number for CO2'
      vname_bec2_diag_2d(3,indxSCHMIDTCO2)=' '
      vname_bec2_diag_2d(4,indxSCHMIDTCO2)=' '

      vname_bec2_diag_2d(1,indxPVO2)='PV_O2'
      vname_bec2_diag_2d(2,indxPVO2)='Piston velocity for O2'
      vname_bec2_diag_2d(3,indxPVO2)='m/s '
      vname_bec2_diag_2d(4,indxPVO2)=' '
      vname_bec2_diag_2d(1,indxPVCO2)='PV_CO2'
      vname_bec2_diag_2d(2,indxPVCO2)='Piston velocity for CO2'
      vname_bec2_diag_2d(3,indxPVCO2)='m/s '
      vname_bec2_diag_2d(4,indxPVCO2)=' '
#ifdef N2O_NEV
      vname_bec2_diag_2d(1,indxfgn2o_nev)='FG_N2O_NEV'
      vname_bec2_diag_2d(2,indxfgn2o_nev)='Air-sea N2O flux (Nevison)'
      vname_bec2_diag_2d(3,indxfgn2o_nev)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o_nev)=' '
#endif
# ifdef Ncycle_SY
# ifdef N2O_TRACER_DECOMP
      vname_bec2_diag_2d(1,indxfgn2o_ao1)='FG_N2O_AO1'
      vname_bec2_diag_2d(2,indxfgn2o_ao1)='Air-sea flux of N2O from AO1'
      vname_bec2_diag_2d(3,indxfgn2o_ao1)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxfgn2o_ao1)=' '

      vname_bec2_diag_2d(1,indxfgn2o_siden)='FG_N2O_SIDEN'
      vname_bec2_diag_2d(2,indxfgn2o_siden)=
     &   'Air-sea flux of N2O from denitrification sink'
      vname_bec2_diag_2d(3,indxfgn2o_siden)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o_siden)=' '

      vname_bec2_diag_2d(1,indxfgn2o_soden)='FG_N2O_SODEN'
      vname_bec2_diag_2d(2,indxfgn2o_soden)=
     &   'Air-sea flux of N2O from denitrification source'
      vname_bec2_diag_2d(3,indxfgn2o_soden)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o_soden)=' '

      vname_bec2_diag_2d(1,indxfgn2o_atm)='FG_N2O_ATM'
      vname_bec2_diag_2d(2,indxfgn2o_atm)=
     &   'Air-sea flux of atmopsheric N2O'
      vname_bec2_diag_2d(3,indxfgn2o_atm)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o_atm)=' '
# endif

      vname_bec2_diag_2d(1,indxfgn2o)='FG_N2O'
      vname_bec2_diag_2d(2,indxfgn2o)=
     &   'Air-sea flux of N2O '
      vname_bec2_diag_2d(3,indxfgn2o)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o)=' '

      vname_bec2_diag_2d(1,indxschmidt_n2o)='SCHMIDT_N2O'
      vname_bec2_diag_2d(2,indxschmidt_n2o)='Schmidt number for N2O'
      vname_bec2_diag_2d(3,indxschmidt_n2o)=' '
      vname_bec2_diag_2d(4,indxschmidt_n2o)=' '

      vname_bec2_diag_2d(1,indxn2osat)='N2OSAT'
      vname_bec2_diag_2d(2,indxn2osat)='N2O saturation concentration'
      vname_bec2_diag_2d(3,indxn2osat)='mmol/m3 '
      vname_bec2_diag_2d(4,indxn2osat)=' '

      vname_bec2_diag_2d(1,indxpvn2o)='PV_N2O'
      vname_bec2_diag_2d(2,indxpvn2o)='Piston velocity for N2O'
      vname_bec2_diag_2d(3,indxpvn2o)='m/s '
      vname_bec2_diag_2d(4,indxpvn2o)=' '

      vname_bec2_diag_2d(1,indxfgn2)='FG_N2'
      vname_bec2_diag_2d(2,indxfgn2)=
     &   'Air-sea flux of N2 '
      vname_bec2_diag_2d(3,indxfgn2)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2)=' '

      vname_bec2_diag_2d(1,indxschmidt_n2)='SCHMIDT_N2'
      vname_bec2_diag_2d(2,indxschmidt_n2)='Schmidt number for N2'
      vname_bec2_diag_2d(3,indxschmidt_n2)=' '
      vname_bec2_diag_2d(4,indxschmidt_n2)=' '

      vname_bec2_diag_2d(1,indxn2sat)='N2SAT'
      vname_bec2_diag_2d(2,indxn2sat)='N2 saturation concentration'
      vname_bec2_diag_2d(3,indxn2sat)='mmol/m3 '
      vname_bec2_diag_2d(4,indxn2sat)=' '

      vname_bec2_diag_2d(1,indxpvn2)='PV_N2'
      vname_bec2_diag_2d(2,indxpvn2)='Piston velocity for N2'
      vname_bec2_diag_2d(3,indxpvn2)='m/s '
      vname_bec2_diag_2d(4,indxpvn2)=' '
#endif
#  if defined CCHEM_MOCSY && defined CCHEM_TODEPTH
      vname_bec2_diag_3d(1,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_3d(2,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_3d(3,indxCO2STAR)='mmol/m3 '
      vname_bec2_diag_3d(4,indxCO2STAR)=' '
#  else
      vname_bec2_diag_2d(1,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_2d(2,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_2d(3,indxCO2STAR)='mmol/m3 '
      vname_bec2_diag_2d(4,indxCO2STAR)=' '
#  endif
      vname_bec2_diag_2d(1,indxDCO2STAR)='DCO2STAR'
      vname_bec2_diag_2d(2,indxDCO2STAR)='DCO2STAR'
      vname_bec2_diag_2d(3,indxDCO2STAR)='mmol/m3'
      vname_bec2_diag_2d(4,indxDCO2STAR)=' '

      vname_bec2_diag_2d(1,indxPCO2OC)='PCO2OC'
      vname_bec2_diag_2d(2,indxPCO2OC)='PCO2OC'
      vname_bec2_diag_2d(3,indxPCO2OC)='not looked up yet '
      vname_bec2_diag_2d(4,indxPCO2OC)=' '

      vname_bec2_diag_2d(1,indxIRONFLUX)='IRON_FLUX'
      vname_bec2_diag_2d(2,indxIRONFLUX)='Iron surface flux'
      vname_bec2_diag_2d(3,indxIRONFLUX)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxIRONFLUX)=' '

      vname_bec2_diag_2d(1,indxSEDDENITRIF)='SED_DENITRIF'
      vname_bec2_diag_2d(2,indxSEDDENITRIF)='Sediment denitrification'
      vname_bec2_diag_2d(3,indxSEDDENITRIF)='mmol/m3/s'
      vname_bec2_diag_2d(4,indxSEDDENITRIF)=' '

#  ifdef CCHEM_MOCSY
#   ifdef CCHEM_TODEPTH
      vname_bec2_diag_3d(1,indxCO3)='CO3'
      vname_bec2_diag_3d(2,indxCO3)='carbonate ion (CO3-) concentration'
      vname_bec2_diag_3d(3,indxCO3)='mmol/m3'
      vname_bec2_diag_3d(4,indxCO3)=' '
      vname_bec2_diag_3d(1,indxHCO3)='HCO3'
      vname_bec2_diag_3d(2,indxHCO3)=
     &  'bicarbonate ion (HCO3-) concentration'
      vname_bec2_diag_3d(3,indxHCO3)='mmol/m3'
      vname_bec2_diag_3d(4,indxHCO3)=' '
      vname_bec2_diag_3d(1,indxOMEGACALC)='OMEGA_CALC'
      vname_bec2_diag_3d(2,indxOMEGACALC)=
     &  'Calcite saturation state'
      vname_bec2_diag_3d(3,indxOMEGACALC)=' '
      vname_bec2_diag_3d(4,indxOMEGACALC)=' '
      vname_bec2_diag_3d(1,indxOMEGAARAG)='OMEGA_ARAG'
      vname_bec2_diag_3d(2,indxOMEGAARAG)=
     &  'Aragonite saturation state'
      vname_bec2_diag_3d(3,indxOMEGAARAG)=' '
      vname_bec2_diag_3d(4,indxOMEGAARAG)=' '
#   else
      vname_bec2_diag_2d(1,indxCO3)='CO3'
      vname_bec2_diag_2d(2,indxCO3)='carbonate ion (CO3-) concentration'
      vname_bec2_diag_2d(3,indxCO3)='mmol/m3'
      vname_bec2_diag_2d(4,indxCO3)=' '
      vname_bec2_diag_2d(1,indxHCO3)='HCO3'
      vname_bec2_diag_2d(2,indxHCO3)=
     &  'bicarbonate ion (HCO3-) concentration'
      vname_bec2_diag_2d(3,indxHCO3)='mmol/m3'
      vname_bec2_diag_2d(4,indxHCO3)=' '
#   endif /* CCHEM_TODEPTH */
#  endif /* CCHEM_MOCSY */

      vname_bec2_diag_2d(1,indxFESEDFLUX)='SED_FE_FLUX'
      vname_bec2_diag_2d(2,indxFESEDFLUX)='Benthic Iron Flux'
      vname_bec2_diag_2d(3,indxFESEDFLUX)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFESEDFLUX)=' '

      vname_bec2_diag_2d(1,indxFLUXTOSED)='FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxFLUXTOSED)='POC reaching the sediments'
      vname_bec2_diag_2d(3,indxFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxCACO3FLUXTOSED)='CACO3_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxCACO3FLUXTOSED)='CaCO3 reaching the sed'
      vname_bec2_diag_2d(3,indxCACO3FLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxCACO3FLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxSIO2FLUXTOSED)='SIO2_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxSIO2FLUXTOSED)='Opal reaching the sed'
      vname_bec2_diag_2d(3,indxSIO2FLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxSIO2FLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxPIRONFLUXTOSED)='PIRON_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxPIRONFLUXTOSED)='P iron reaching the sed'
      vname_bec2_diag_2d(3,indxPIRONFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxPIRONFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxDUSTFLUXTOSED)='DUST_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxDUSTFLUXTOSED)='Dust reaching the sed'
      vname_bec2_diag_2d(3,indxDUSTFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxDUSTFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxPOCSEDLOSS)='POC_SED_LOSS'
      vname_bec2_diag_2d(2,indxPOCSEDLOSS)='POC lost to sed at bottom'
      vname_bec2_diag_2d(3,indxPOCSEDLOSS)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxPOCSEDLOSS)=' '

      vname_bec2_diag_2d(1,indxOTHERREMIN)='OTHER_REMIN'
      vname_bec2_diag_2d(2,indxOTHERREMIN)='other remin pathways in sed'
      vname_bec2_diag_2d(3,indxOTHERREMIN)='mmol/m3/s'
      vname_bec2_diag_2d(4,indxOTHERREMIN)=' '

      vname_bec2_diag_2d(1,indxCACO3SEDLOSS)='P_CACO3_SED_LOSS'
      vname_bec2_diag_2d(2,indxCACO3SEDLOSS)='CaCO3 lost to sed'
      vname_bec2_diag_2d(3,indxCACO3SEDLOSS)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxCACO3SEDLOSS)=' '

      vname_bec2_diag_2d(1,indxSIO2SEDLOSS)='SIO2_SED_LOSS'
      vname_bec2_diag_2d(2,indxSIO2SEDLOSS)='SiO2 lost to sed'
      vname_bec2_diag_2d(3,indxSIO2SEDLOSS)='mmol Si/m2/s'
      vname_bec2_diag_2d(4,indxSIO2SEDLOSS)=' '


      ! MF: Control BEC2_DIAG Output vars
      wrtavg_bec2_diag_2d=(/.true.,.false.,.false.,.false.,.false.,   ! pCO2air, PARinc,FG_O2,FG_CO2, WS10m
     &                     .false.,.false.,.false.,.false.,.false.,   ! XKW, ATM_PRESS,SCHMIDT_O2,O2SAT, SCHMIDT_CO2
     &                     .false.,.false.,.true.,.true.,.false.,    ! PVO2, PV_CO2,IRON_FLUX, SED_DENITRI, PH
     &                     .true.,.false.,.false.,.false.,          ! PCO2, CO2STAR, PCO2OC, DCO2STAR
#ifdef Ncycle_SY
     &                     .false.,.false.,.true.,.true., !SCHMIDT_N2O, PV_N2O, N2OSAT, FG_N2O
     &                     .false.,.false.,.true.,.true., !SCHMIDT_N2,PV_N2, N2SAT, FG_N2
# ifdef N2O_TRACER_DECOMP
     &                     .true.,.true.,.true.,.true., !FG_N2O_AO1,  FG_N2O_SIDEN, FG_N2O_SODEN, FG_N2O_ATM
# endif
# endif
# ifdef N2O_NEV
     &                     .true.,!FG_N2O_NEV
#endif
     &                     .false.,.false.,.false.,.false.,.false.,   ! FE_SED_FLUX,FLUX_TO_SED,CACO3_FLUX_TO_SED,SIO2_FLUX_TO_SED,PIRON_FLUX_TO_SED
     &                     .false.,.true.,.true.,.false.,.false./)  ! DUST_FLUX_SED, POC_SED_LOSS, OTHER_REMIN,CaCO3_sed_loss,P_SiO2_sed_loss

      !
      ! More 3d diagnostics:
      !
      vname_bec2_diag_3d(1,indxPAR)='PAR'
      vname_bec2_diag_3d(2,indxPAR)=
     &  'Photosynthetically active radiation'
      vname_bec2_diag_3d(3,indxPAR)='W m-2'
      vname_bec2_diag_3d(4,indxPAR)='  '

      vname_bec2_diag_3d(1,indxPOCFLUXIN)='POC_FLUX_IN'
      vname_bec2_diag_3d(2,indxPOCFLUXIN)='POC flux into cell'
      vname_bec2_diag_3d(3,indxPOCFLUXIN)='mmol C/m2/s '
      vname_bec2_diag_3d(4,indxPOCFLUXIN)=' '

      vname_bec2_diag_3d(1,indxPOCPROD)='POC_PROD'
      vname_bec2_diag_3d(2,indxPOCPROD)='POC production'
      vname_bec2_diag_3d(3,indxPOCPROD)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPOCPROD)=' '
      vname_bec2_diag_3d(1,indxPOCREMIN)='POC_REMIN'
      vname_bec2_diag_3d(2,indxPOCREMIN)='POC remineralization'
      vname_bec2_diag_3d(3,indxPOCREMIN)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxPOCREMIN)=' '

      vname_bec2_diag_3d(1,indxCACO3FLUXIN)='CACO3_FLUX_IN'
      vname_bec2_diag_3d(2,indxCACO3FLUXIN)=
     &  'Incoming flux of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxCACO3FLUXIN)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,indxCACO3FLUXIN)=' '
      vname_bec2_diag_3d(1,indxPCACO3PROD)='PCACO3_PROD'
      vname_bec2_diag_3d(2,indxPCACO3PROD)=
     &  'Production of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxPCACO3PROD)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxPCACO3PROD)=' '
      vname_bec2_diag_3d(1,indxCACO3REMIN)='CACO3_REMIN'
      vname_bec2_diag_3d(2,indxCACO3REMIN)=
     &  'Remineralization of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxCACO3REMIN)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3REMIN)=' '

      vname_bec2_diag_3d(1,indxSIO2FLUXIN)='SIO2_FLUX_IN'
      vname_bec2_diag_3d(2,indxSIO2FLUXIN)=
     &  'Incoming flux of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2FLUXIN)='mmol Si/m2/s '
      vname_bec2_diag_3d(4,indxSIO2FLUXIN)=' '
      vname_bec2_diag_3d(1,indxSIO2PROD)='SIO2_PROD'
      vname_bec2_diag_3d(2,indxSIO2PROD)=
     &  'Production of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2PROD)='mmol Si/m3/s '
      vname_bec2_diag_3d(4,indxSIO2PROD)=' '
      vname_bec2_diag_3d(1,indxSIO2REMIN)='SIO2_REMIN'
      vname_bec2_diag_3d(2,indxSIO2REMIN)=
     &  'Remineralization of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2REMIN)='mmol Si/m3/s'
      vname_bec2_diag_3d(4,indxSIO2REMIN)=' '

      vname_bec2_diag_3d(1,indxDUSTFLUXIN)='DUST_FLUX_IN'
      vname_bec2_diag_3d(2,indxDUSTFLUXIN)=
     &  'Incoming dust flux'
      vname_bec2_diag_3d(3,indxDUSTFLUXIN)='kg dust/m2/s '
      vname_bec2_diag_3d(4,indxDUSTFLUXIN)=' '
      vname_bec2_diag_3d(1,indxDUSTREMIN)='DUST_REMIN'
      vname_bec2_diag_3d(2,indxDUSTREMIN)=
     &  'Remineralization of dust'
      vname_bec2_diag_3d(3,indxDUSTREMIN)='kg dust/m3/s '
      vname_bec2_diag_3d(4,indxDUSTREMIN)=' '

      vname_bec2_diag_3d(1,indxPIRONFLUXIN)='P_IRON_FLUX_IN'
      vname_bec2_diag_3d(2,indxPIRONFLUXIN)=
     &  'Incoming flux of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONFLUXIN)='mmol Fe/m2/s '
      vname_bec2_diag_3d(4,indxPIRONFLUXIN)=' '
      vname_bec2_diag_3d(1,indxPIRONPROD)='P_IRON_PROD'
      vname_bec2_diag_3d(2,indxPIRONPROD)=
     &  'Production of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONPROD)='mmol Fe/m3/s '
      vname_bec2_diag_3d(4,indxPIRONPROD)=' '
      vname_bec2_diag_3d(1,indxPIRONREMIN)='P_IRON_REMIN'
      vname_bec2_diag_3d(2,indxPIRONREMIN)=
     &  'Remineralization of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONREMIN)='mmol Fe/m3/s '
      vname_bec2_diag_3d(4,indxPIRONREMIN)=' '

      vname_bec2_diag_3d(1,indxGRAZESP)='GRAZE_SP'
      vname_bec2_diag_3d(2,indxGRAZESP)=
     &  'Grazing rate on small phytoplankton'
      vname_bec2_diag_3d(3,indxGRAZESP)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZESP)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAT)='GRAZE_DIAT'
      vname_bec2_diag_3d(2,indxGRAZEDIAT)=
     &  'Grazing rate on diatoms'
      vname_bec2_diag_3d(3,indxGRAZEDIAT)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAT)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAZ)='GRAZE_DIAZ'
      vname_bec2_diag_3d(2,indxGRAZEDIAZ)=
     &  'Grazing rate on diazotrophs'
      vname_bec2_diag_3d(3,indxGRAZEDIAZ)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAZ)=' '

      vname_bec2_diag_3d(1,indxSPLOSS)='SP_LOSS'
      vname_bec2_diag_3d(2,indxSPLOSS)=
     &  'Small phytoplankton non-grazing mortality'
      vname_bec2_diag_3d(3,indxSPLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxSPLOSS)=' '
      vname_bec2_diag_3d(1,indxDIATLOSS)='DIAT_LOSS'
      vname_bec2_diag_3d(2,indxDIATLOSS)=
     &  'Diatom non-grazing mortality'
      vname_bec2_diag_3d(3,indxDIATLOSS)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxDIATLOSS)=' '
      vname_bec2_diag_3d(1,indxZOOLOSS)='ZOO_LOSS'
      vname_bec2_diag_3d(2,indxZOOLOSS)=
     &  'Mortality due to higher trophic grazing on zooplankton'
      vname_bec2_diag_3d(3,indxZOOLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxZOOLOSS)=' '

      vname_bec2_diag_3d(1,indxSPAGG)='SP_AGG'
      vname_bec2_diag_3d(2,indxSPAGG)=
     &  'Aggregation of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxSPAGG)=' '
      vname_bec2_diag_3d(1,indxDIATAGG)='DIAT_AGG'
      vname_bec2_diag_3d(2,indxDIATAGG)=
     &  'Aggregation of diatoms'
      vname_bec2_diag_3d(3,indxDIATAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDIATAGG)=' '

      vname_bec2_diag_3d(1,indxPHOTOCSP)='PHOTOC_SP'
      vname_bec2_diag_3d(2,indxPHOTOCSP)=
     &  'C fixation rate by small phytoplankton'
      vname_bec2_diag_3d(3,indxPHOTOCSP)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCSP)=' '
      vname_bec2_diag_3d(1,indxPHOTOCDIAT)='PHOTOC_DIAT'
      vname_bec2_diag_3d(2,indxPHOTOCDIAT)=
     &  'C fixation rate by diatoms'
      vname_bec2_diag_3d(3,indxPHOTOCDIAT)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCDIAT)=' '

      vname_bec2_diag_3d(1,indxTOTPROD)='TOT_PROD'
      vname_bec2_diag_3d(2,indxTOTPROD)=
     &  'Total autotroph production'
      vname_bec2_diag_3d(3,indxTOTPROD)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxTOTPROD)=' '
      vname_bec2_diag_3d(1,indxDOCPROD)='DOC_PROD'
      vname_bec2_diag_3d(2,indxDOCPROD)='DOC production'
      vname_bec2_diag_3d(3,indxDOCPROD)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDOCPROD)=' '

      vname_bec2_diag_3d(1,indxDOCREMIN)='DOC_REMIN'
      vname_bec2_diag_3d(2,indxDOCREMIN)=
     &  'Remineralization rate of DOC'
      vname_bec2_diag_3d(3,indxDOCREMIN)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxDOCREMIN)=' '
      vname_bec2_diag_3d(1,indxFESCAVENGE)='FE_SCAVENGE'
      vname_bec2_diag_3d(2,indxFESCAVENGE)=
     &  'Loss of dissolved iron, scavenging'
      vname_bec2_diag_3d(3,indxFESCAVENGE)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxFESCAVENGE)=' '

      vname_bec2_diag_3d(1,indxSPNLIM)='SP_N_LIM'
      vname_bec2_diag_3d(2,indxSPNLIM)=
     &  'N limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPNLIM)=' '
      vname_bec2_diag_3d(4,indxSPNLIM)=' '
      vname_bec2_diag_3d(1,indxSPFEUPTAKE)='SP_FE_LIM'
      vname_bec2_diag_3d(2,indxSPFEUPTAKE)=
     &  'Fe limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxSPFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxSPPO4UPTAKE)='SP_PO4_LIM'
      vname_bec2_diag_3d(2,indxSPPO4UPTAKE)=
     &  'PO4 limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxSPPO4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxSPLIGHTLIM)='SP_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxSPLIGHTLIM)=
     &  'Small phytoplankton light limitation factor'
      vname_bec2_diag_3d(3,indxSPLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxSPLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxDIATNLIM)='DIAT_N_LIM'
      vname_bec2_diag_3d(2,indxDIATNLIM)=
     &  'N limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATNLIM)=' '
      vname_bec2_diag_3d(4,indxDIATNLIM)=' '
      vname_bec2_diag_3d(1,indxDIATFEUPTAKE)='DIAT_FE_LIM'
      vname_bec2_diag_3d(2,indxDIATFEUPTAKE)=
     &  'Fe limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATPO4UPTAKE)='DIAT_PO4_LIM'
      vname_bec2_diag_3d(2,indxDIATPO4UPTAKE)=
     &  'PO4 limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATPO4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATSIO3UPTAKE)='DIAT_SIO3_LIM'
      vname_bec2_diag_3d(2,indxDIATSIO3UPTAKE)=
     &  'SiO3 limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATSIO3UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATSIO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIATLIGHTLIM)='DIAT_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxDIATLIGHTLIM)=
     &  'Diatom light limitation factor'
      vname_bec2_diag_3d(3,indxDIATLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxDIATLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxCACO3PROD)='CACO3_PROD'
      vname_bec2_diag_3d(2,indxCACO3PROD)=
     &  'Production of CaCO3 by small phytoplankton'
      vname_bec2_diag_3d(3,indxCACO3PROD)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3PROD)=' '
      vname_bec2_diag_3d(1,indxDIAZNFIX)='DIAZ_NFIX'
      vname_bec2_diag_3d(2,indxDIAZNFIX)=
     &  'Total nitrogen fixation by diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNFIX)='mmol/m3/s   '
      vname_bec2_diag_3d(4,indxDIAZNFIX)=' '

      vname_bec2_diag_3d(1,indxDIAZLOSS)='DIAZ_LOSS'
      vname_bec2_diag_3d(2,indxDIAZLOSS)=
     &  'Diazotroph non-grazing mortality'
      vname_bec2_diag_3d(3,indxDIAZLOSS)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZLOSS)=' '
      vname_bec2_diag_3d(1,indxPHOTOCDIAZ)='PHOTOC_DIAZ'
      vname_bec2_diag_3d(2,indxPHOTOCDIAZ)=
     &   'Diazotroph C-fixation'
      vname_bec2_diag_3d(3,indxPHOTOCDIAZ)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxPHOTOCDIAZ)=' '

      vname_bec2_diag_3d(1,indxDIAZPO4UPTAKE)='DIAZ_PO4_LIM'
      vname_bec2_diag_3d(2,indxDIAZPO4UPTAKE)=
     &  'PO4 limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIAZPO4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZFEUPTAKE)='DIAZ_FE_LIM'
      vname_bec2_diag_3d(2,indxDIAZFEUPTAKE)=
     &  'Fe limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIAZFEUPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIAZLIGHTLIM)='DIAZ_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxDIAZLIGHTLIM)=
     &   'Diazotroph light limitation'
      vname_bec2_diag_3d(3,indxDIAZLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxDIAZLIGHTLIM)=' '
      vname_bec2_diag_3d(1,indxFESCAVENGERATE)='FE_SCAVENGE_RATE'
      vname_bec2_diag_3d(2,indxFESCAVENGERATE)=
     &  'Annual scavenging rate of iron as % of ambient'
      vname_bec2_diag_3d(3,indxFESCAVENGERATE)=' '
      vname_bec2_diag_3d(4,indxFESCAVENGERATE)=' '

      vname_bec2_diag_3d(1,indxDONPROD)='DON_PROD'
      vname_bec2_diag_3d(2,indxDONPROD)=
     &  'Production of dissolved organic N'
      vname_bec2_diag_3d(3,indxDONPROD)=' '
      vname_bec2_diag_3d(4,indxDONPROD)=' '
      vname_bec2_diag_3d(1,indxDONREMIN)='DON_REMIN'
      vname_bec2_diag_3d(2,indxDONREMIN)=
     &  'Remineralization rate of DON'
      vname_bec2_diag_3d(3,indxDONREMIN)='mmol N/m3/s  '
      vname_bec2_diag_3d(4,indxDONREMIN)=' '

      vname_bec2_diag_3d(1,indxDOFEPROD)='DOFE_PROD'
      vname_bec2_diag_3d(2,indxDOFEPROD)=
     &  'Production of dissolved organic iron'
      vname_bec2_diag_3d(3,indxDOFEPROD)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOFEPROD)=' '
      vname_bec2_diag_3d(1,indxDOFEREMIN)='DOFE_REMIN'
      vname_bec2_diag_3d(2,indxDOFEREMIN)=
     &  'Remineralization rate of DOFE'
      vname_bec2_diag_3d(3,indxDOFEREMIN)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOFEREMIN)=' '

      vname_bec2_diag_3d(1,indxDOPPROD)='DOP_PROD'
      vname_bec2_diag_3d(2,indxDOPPROD)=
     &  'Production of dissolved organic P'
      vname_bec2_diag_3d(3,indxDOPPROD)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOPPROD)=' '
      vname_bec2_diag_3d(1,indxDOPREMIN)='DOP_REMIN'
      vname_bec2_diag_3d(2,indxDOPREMIN)=
     &  'Remineralization rate of DOP'
      vname_bec2_diag_3d(3,indxDOPREMIN)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxDOPREMIN)=' '

      vname_bec2_diag_3d(1,indxDIATSIUPTAKE)='DIAT_SI_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATSIUPTAKE)=
     &  'Silicon uptake rate by diatoms'
      vname_bec2_diag_3d(3,indxDIATSIUPTAKE)='mmol Si/m3/s  '
      vname_bec2_diag_3d(4,indxDIATSIUPTAKE)=' '

      vname_bec2_diag_3d(1,indxIRONUPTAKESP)='PHOTOFE_SP'
      vname_bec2_diag_3d(2,indxIRONUPTAKESP)=
     &  'Iron uptake rate by small phytoplankton'
      vname_bec2_diag_3d(3,indxIRONUPTAKESP)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKESP)=' '
      vname_bec2_diag_3d(1,indxIRONUPTAKEDIAT)='PHOTOFE_DIAT'
      vname_bec2_diag_3d(2,indxIRONUPTAKEDIAT)=
     &  'Iron uptake rate by diatoms'
      vname_bec2_diag_3d(3,indxIRONUPTAKEDIAT)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKEDIAT)=' '
      vname_bec2_diag_3d(1,indxIRONUPTAKEDIAZ)='PHOTOFE_DIAZ'
      vname_bec2_diag_3d(2,indxIRONUPTAKEDIAZ)=
     &  'Iron uptake rate by diatotrophs'
      vname_bec2_diag_3d(3,indxIRONUPTAKEDIAZ)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKEDIAZ)=' '

      vname_bec2_diag_3d(1,indxNITRIF)='NITRIF'
      vname_bec2_diag_3d(2,indxNITRIF)='Nitrification'
      vname_bec2_diag_3d(3,indxNITRIF)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxNITRIF)=' '
      vname_bec2_diag_3d(1,indxDENITRIF)='DENITRIF'
      vname_bec2_diag_3d(2,indxDENITRIF)='Denitrification'
      vname_bec2_diag_3d(3,indxDENITRIF)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDENITRIF)=' '

      vname_bec2_diag_3d(1,indxSPNO3UPTAKE)='SP_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNO3UPTAKE)=
     &   'NO3 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNO3UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATNO3UPTAKE)='DIAT_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNO3UPTAKE)=
     &   'NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNO3UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZNO3UPTAKE)='DIAZ_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNO3UPTAKE)=
     &   'NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxSPNH4UPTAKE)='SP_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNH4UPTAKE)=
     &   'NH4 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNH4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATNH4UPTAKE)='DIAT_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNH4UPTAKE)=
     &   'NH4 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNH4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZNH4UPTAKE)='DIAZ_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNH4UPTAKE)=
     &   'NH4 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNH4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxGRAZEDICSP)='SP_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICSP)=
     &   'Small phyto grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICSP)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICSP)=' '
      vname_bec2_diag_3d(1,indxGRAZEDICDIAT)='DIAT_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICDIAT)=
     &   'Diatom grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICDIAT)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICDIAT)=' '
      vname_bec2_diag_3d(1,indxGRAZEDICDIAZ)='DIAZ_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICDIAZ)=
     &   'Diazotroph grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICDIAZ)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICDIAZ)=' '

      vname_bec2_diag_3d(1,indxLOSSDICSP)='SP_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICSP)=
     &   'Small phyto non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICSP)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICSP)=' '
      vname_bec2_diag_3d(1,indxLOSSDICDIAT)='DIAT_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICDIAT)=
     &   'Diatom non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICDIAT)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICDIAT)=' '
      vname_bec2_diag_3d(1,indxLOSSDICDIAZ)='DIAZ_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICDIAZ)=
     &   'Diazotroph non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICDIAZ)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICDIAZ)=' '

      vname_bec2_diag_3d(1,indxZOOLOSSDIC)='ZOO_LOSS_DIC'
      vname_bec2_diag_3d(2,indxZOOLOSSDIC)=
     &   'Zooplankton loss routed to DIC'
      vname_bec2_diag_3d(3,indxZOOLOSSDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxZOOLOSSDIC)=' '

      vname_bec2_diag_3d(1,indxDIAZAGG)='DIAZ_AGG'
      vname_bec2_diag_3d(2,indxDIAZAGG)=
     &  'Aggregation of diatoms'
      vname_bec2_diag_3d(3,indxDIAZAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDIAZAGG)=' '

      vname_bec2_diag_3d(1,indxGRAZESPZOO)='GRAZE_SP_ZOO'
      vname_bec2_diag_3d(2,indxGRAZESPZOO)=
     &  'Grazing rate on small phytoplankton routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZESPZOO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZESPZOO)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIATZOO)='GRAZE_DIAT_ZOO'
      vname_bec2_diag_3d(2,indxGRAZEDIATZOO)=
     &  'Grazing rate on diatoms routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZEDIATZOO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIATZOO)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAZZOO)='GRAZE_DIAZ_ZOO'
      vname_bec2_diag_3d(2,indxGRAZEDIAZZOO)=
     &  'Grazing rate on diazotrophs routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZEDIAZZOO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAZZOO)=' '

      vname_bec2_diag_3d(1,indxSPQCACO3)='QCACO3'
      vname_bec2_diag_3d(2,indxSPQCACO3)=
     &  'Small phyto CaCO3/C ratio'
      vname_bec2_diag_3d(3,indxSPQCACO3)='mmol CaCO3/mmol C'
      vname_bec2_diag_3d(4,indxSPQCACO3)=' '

      vname_bec2_diag_3d(1,indxSPPHOTOACC)='PHOTOACC_SP'
      vname_bec2_diag_3d(2,indxSPPHOTOACC)=
     &  'Chl synthesis in photoadaptation for small phyto'
      vname_bec2_diag_3d(3,indxSPPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxSPPHOTOACC)=' '
      vname_bec2_diag_3d(1,indxDIATPHOTOACC)='PHOTOACC_DIAT'
      vname_bec2_diag_3d(2,indxDIATPHOTOACC)=
     &  'Chl synthesis in photoadaptation for diatoms'
      vname_bec2_diag_3d(3,indxDIATPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxDIATPHOTOACC)=' '
      vname_bec2_diag_3d(1,indxDIAZPHOTOACC)='PHOTOACC_DIAZ'
      vname_bec2_diag_3d(2,indxDIAZPHOTOACC)=
     &  'Chl synthesis in photoadaptation for diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxDIAZPHOTOACC)=' '

      vname_bec2_diag_3d(1,indxSPCZERO)='SPC_ZERO'
      vname_bec2_diag_3d(2,indxSPCZERO)=
     &  'Change caused by setting negative SPC values to zero'
      vname_bec2_diag_3d(3,indxSPCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxSPCZERO)=' '
      vname_bec2_diag_3d(1,indxDIATCZERO)='DIATC_ZERO'
      vname_bec2_diag_3d(2,indxDIATCZERO)=
     &  'Change caused by setting negative DIATC values to zero'
      vname_bec2_diag_3d(3,indxDIATCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDIATCZERO)=' '
      vname_bec2_diag_3d(1,indxDIAZCZERO)='DIAZC_ZERO'
      vname_bec2_diag_3d(2,indxDIAZCZERO)=
     &  'Change caused by setting negative DIAZC values to zero'
      vname_bec2_diag_3d(3,indxDIAZCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDIAZCZERO)=' '
      vname_bec2_diag_3d(1,indxDOCZERO)='DOC_ZERO'
      vname_bec2_diag_3d(2,indxDOCZERO)=
     &  'Change caused by setting negative DOC values to zero'
      vname_bec2_diag_3d(3,indxDOCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDOCZERO)=' '
      vname_bec2_diag_3d(1,indxZOOCZERO)='ZOOC_ZERO'
      vname_bec2_diag_3d(2,indxZOOCZERO)=
     &  'Change caused by setting negative ZOOC values to zero'
      vname_bec2_diag_3d(3,indxZOOCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxZOOCZERO)=' '
      vname_bec2_diag_3d(1,indxSPCACO3ZERO)='SPCACO3_ZERO'
      vname_bec2_diag_3d(2,indxSPCACO3ZERO)=
     &  'Change caused by setting negative SPCACO3 values to zero'
      vname_bec2_diag_3d(3,indxSPCACO3ZERO)='mMol CaCO3/m3/s'
      vname_bec2_diag_3d(4,indxSPCACO3ZERO)=' '

      vname_bec2_diag_3d(1,indxDONRREMIN) = 'DONr_REMIN'
      vname_bec2_diag_3d(2,indxDONRREMIN)=
     &  'Portion of refractory DON remineralized'
      vname_bec2_diag_3d(3,indxDONRREMIN)='mMol N/m3/s'
      vname_bec2_diag_3d(4,indxDONRREMIN)=' '

      vname_bec2_diag_3d(1,indxTOTCHL)='TOT_CHL'
      vname_bec2_diag_3d(2,indxTOTCHL)=
     &  'Total Chlorophyll'
      vname_bec2_diag_3d(3,indxTOTCHL)='mg Chl/m3'
      vname_bec2_diag_3d(4,indxTOTCHL)=' '

      vname_bec2_diag_3d(1,indxTOTPHYTOC)='TOT_PHYTOC'
      vname_bec2_diag_3d(2,indxTOTPHYTOC)=
     &  'Total Phytoplankton Carbon'
      vname_bec2_diag_3d(3,indxTOTPHYTOC)='mMol C/m3'
      vname_bec2_diag_3d(4,indxTOTPHYTOC)=' '

      vname_bec2_diag_3d(1,indxSPPLIM)='SP_P_LIM'
      vname_bec2_diag_3d(2,indxSPPLIM)=
     &  'Small phytoplankton P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxSPPLIM)=' '
      vname_bec2_diag_3d(4,indxSPPLIM)=' '

      vname_bec2_diag_3d(1,indxDIATPLIM)='DIAT_P_LIM'
      vname_bec2_diag_3d(2,indxDIATPLIM)=
     &  'Diatom P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxDIATPLIM)=' '
      vname_bec2_diag_3d(4,indxDIATPLIM)=' '

      vname_bec2_diag_3d(1,indxDIAZPLIM)='DIAZ_P_LIM'
      vname_bec2_diag_3d(2,indxDIAZPLIM)=
     &  'Diazotroph P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxDIAZPLIM)=' '
      vname_bec2_diag_3d(4,indxDIAZPLIM)=' '

      vname_bec2_diag_3d(1,indxo2cons)='O2_CONSUMPTION'
      vname_bec2_diag_3d(2,indxo2cons)=
     &  'O2 consumption rate'
      vname_bec2_diag_3d(3,indxo2cons)='mmol O2/m3/s'
      vname_bec2_diag_3d(4,indxo2cons)=' '

      vname_bec2_diag_3d(1,indxo2prod)='O2_PRODUCTION'
      vname_bec2_diag_3d(2,indxo2prod)=
     &  'O2 production rate'
      vname_bec2_diag_3d(3,indxo2prod)='mmol O2/m3/s'
      vname_bec2_diag_3d(4,indxo2prod)=' '
# ifdef N2O_NEV
      vname_bec2_diag_3d(1,indxn2oconsnev)='N2O_CONS_NEV'
      vname_bec2_diag_3d(2,indxn2oconsnev)=
     &  'N2O consumption rate (Nevison)'
      vname_bec2_diag_3d(3,indxn2oconsnev)='mmol N2O/m3/s'
      vname_bec2_diag_3d(4,indxn2oconsnev)=' '

      vname_bec2_diag_3d(1,indxn2oprodnev)='N2O_PROD_NEV'
      vname_bec2_diag_3d(2,indxn2oprodnev)=
     &  'N2O production rate (Nevison)'
      vname_bec2_diag_3d(3,indxn2oprodnev)='mmol N2O/m3/s'
      vname_bec2_diag_3d(4,indxn2oprodnev)=' '
# endif


# ifdef BEC_COCCO
      vname_bec2_diag_3d(1,indxGRAZECOCCO)='GRAZE_COCCO'
      vname_bec2_diag_3d(2,indxGRAZECOCCO)=
     &  'Grazing rate on Coccos'
      vname_bec2_diag_3d(3,indxGRAZECOCCO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZECOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCOLOSS)='COCCO_LOSS'
      vname_bec2_diag_3d(2,indxCOCCOLOSS)=
     &  'Cocco non-grazing mortality'
      vname_bec2_diag_3d(3,indxCOCCOLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxCOCCOLOSS)=' '

      vname_bec2_diag_3d(1,indxCOCCOAGG)='COCCO_AGG'
      vname_bec2_diag_3d(2,indxCOCCOAGG)=
     &  'Aggregation of coccos'
      vname_bec2_diag_3d(3,indxCOCCOAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxCOCCOAGG)=' '

      vname_bec2_diag_3d(1,indxPHOTOCCOCCO)='PHOTOC_COCCO'
      vname_bec2_diag_3d(2,indxPHOTOCCOCCO)=
     &  'C fixation rate by Coccos'
      vname_bec2_diag_3d(3,indxPHOTOCCOCCO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCCOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCONLIM)='COCCO_N_LIM'
      vname_bec2_diag_3d(2,indxCOCCONLIM)=
     &  'N limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCONLIM)=' '
      vname_bec2_diag_3d(1,indxCOCCOFEUPTAKE)='COCCO_FE_LIM'
      vname_bec2_diag_3d(2,indxCOCCOFEUPTAKE)=
     &  'Fe limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCOFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxCOCCOFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxCOCCOPO4UPTAKE)='COCCO_PO4_LIM'
      vname_bec2_diag_3d(2,indxCOCCOPO4UPTAKE)=
     &  'PO4 limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCOPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxCOCCOPO4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCOLIGHTLIM)='COCCO_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxCOCCOLIGHTLIM)=
     &  'Cocco light limitation factor'
      vname_bec2_diag_3d(3,indxCOCCOLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCOLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxCACO3PRODCOCCO)='CACO3_PROD_COCCO'
      vname_bec2_diag_3d(2,indxCACO3PRODCOCCO)=
     &  'Production of CaCO3 by Coccos'
      vname_bec2_diag_3d(3,indxCACO3PRODCOCCO)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3PRODCOCCO)=' '

      vname_bec2_diag_3d(1,indxIRONUPTAKECOCCO)='PHOTOFE_COCCO'
      vname_bec2_diag_3d(2,indxIRONUPTAKECOCCO)=
     &  'Iron uptake rate by Coccos'
      vname_bec2_diag_3d(3,indxIRONUPTAKECOCCO)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKECOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCONO3UPTAKE)='COCCO_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxCOCCONO3UPTAKE)=
     &   'NO3 uptake rate of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCONO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCONH4UPTAKE)='COCCO_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxCOCCONH4UPTAKE)=
     &   'NH4 uptake rate of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCONH4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCOGRAZEDIC)='COCCO_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxCOCCOGRAZEDIC)=
     &   'Cocco grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxCOCCOGRAZEDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCOGRAZEDIC)=' '

      vname_bec2_diag_3d(1,indxCOCCOLOSSDIC)='COCCO_LOSS_DIC'
      vname_bec2_diag_3d(2,indxCOCCOLOSSDIC)=
     &   'Cocco non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxCOCCOLOSSDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCOLOSSDIC)=' '

      vname_bec2_diag_3d(1,indxGRAZECOCCOZOO)='GRAZE_COCCO_ZOO'
      vname_bec2_diag_3d(2,indxGRAZECOCCOZOO)=
     &  'Grazing rate on Coccos routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZECOCCOZOO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZECOCCOZOO)=' '

      vname_bec2_diag_3d(1,indxQCACO3COCCO)='QCACO3_COCCO'
      vname_bec2_diag_3d(2,indxQCACO3COCCO)=
     &  'Cocco CaCO3/C ratio'
      vname_bec2_diag_3d(3,indxQCACO3COCCO)='mmol CaCO3/mmol C'
      vname_bec2_diag_3d(4,indxQCACO3COCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCOPHOTOACC)='PHOTOACC_COCCO'
      vname_bec2_diag_3d(2,indxCOCCOPHOTOACC)=
     &  'Chl synthesis in photoadaptation for Coccos'
      vname_bec2_diag_3d(3,indxCOCCOPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxCOCCOPHOTOACC)=' '

      vname_bec2_diag_3d(1,indxCOCCOPLIM)='COCCO_P_LIM'
      vname_bec2_diag_3d(2,indxCOCCOPLIM)=
     &  'Cocco P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxCOCCOPLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCOPLIM)=' '

# endif /* BEC_COCCO */
# ifdef Ncycle_SY
      vname_bec2_diag_3d(1,indxammox)='AMMOX'
      vname_bec2_diag_3d(2,indxammox)=
     &  'Ammonium oxidation (NH4 -> NO2)'
      vname_bec2_diag_3d(3,indxammox)='mmol N/m3/s '
      vname_bec2_diag_3d(4,indxammox)=' '

      vname_bec2_diag_3d(1,indxnitrox)='NITROX'
      vname_bec2_diag_3d(2,indxnitrox)=
     &  'Nitrite oxidation (NO2 -> NO3)'
      vname_bec2_diag_3d(3,indxnitrox)='mmol N/m3/s'
      vname_bec2_diag_3d(4,indxnitrox)=' '

      vname_bec2_diag_3d(1,indxanammox)='ANAMMOX'
      vname_bec2_diag_3d(2,indxanammox)=
     &  'Anammox (NO2+NH4 -> N2)'
      vname_bec2_diag_3d(3,indxanammox)='mmol N/m3/s'
      vname_bec2_diag_3d(4,indxanammox)=' '

      vname_bec2_diag_3d(1,indxDENITRIF1)='DENITRIF1'
      vname_bec2_diag_3d(2,indxDENITRIF1)='Denitrification1'
      vname_bec2_diag_3d(3,indxDENITRIF1)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF1)=' '

      vname_bec2_diag_3d(1,indxDENITRIF2)='DENITRIF2'
      vname_bec2_diag_3d(2,indxDENITRIF2)='Denitrification2'
      vname_bec2_diag_3d(3,indxDENITRIF2)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF2)=' '

      vname_bec2_diag_3d(1,indxDENITRIF3)='DENITRIF3'
      vname_bec2_diag_3d(2,indxDENITRIF3)='Denitrification3'
      vname_bec2_diag_3d(3,indxDENITRIF3)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF3)=' '

      vname_bec2_diag_3d(1,indxN2OAMMOX)='N2OAMMOX'
      vname_bec2_diag_3d(2,indxN2OAMMOX)='N2O produced by Ammox'
      vname_bec2_diag_3d(3,indxN2OAMMOX)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OAMMOX)=''

      vname_bec2_diag_3d(1,indxN2OSODEN_CONS)='N2OSODEN_CONS'
      vname_bec2_diag_3d(2,indxN2OSODEN_CONS)='N2O_SODEN consumption'
      vname_bec2_diag_3d(3,indxN2OSODEN_CONS)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OSODEN_CONS)=''

      vname_bec2_diag_3d(1,indxN2OAO1_CONS)='N2OAO1_CONS'
      vname_bec2_diag_3d(2,indxN2OAO1_CONS)='N2O_AO1 consumption'
      vname_bec2_diag_3d(3,indxN2OAO1_CONS)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OAO1_CONS)=''

      vname_bec2_diag_3d(1,indxN2OSIDEN_CONS)='N2OSIDEN_CONS'
      vname_bec2_diag_3d(2,indxN2OSIDEN_CONS)='N2O_SIDEN consumption'
      vname_bec2_diag_3d(3,indxN2OSIDEN_CONS)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OSIDEN_CONS)=''

      vname_bec2_diag_3d(1,indxN2OATM_CONS)='N2OATM_CONS'
      vname_bec2_diag_3d(2,indxN2OATM_CONS)='N2O_ATM consumption'
      vname_bec2_diag_3d(3,indxN2OATM_CONS)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OATM_CONS)=''

      vname_bec2_diag_3d(1,indxSPNO2UPTAKE)='SP_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNO2UPTAKE)=
     &   'NO2 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNO2UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIATNO2UPTAKE)='DIAT_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNO2UPTAKE)=
     &   'NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNO2UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIAZNO2UPTAKE)='DIAZ_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNO2UPTAKE)=
     &   'NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNO2UPTAKE)=' '

# endif
# endif /* BEC2_DIAG */

# ifdef USE_EXPLICIT_VSINK
! For output files:
      vname(1,indxdusthard)='dust_hard                             '
      vname(2,indxdusthard)='Hard component of dust                '
      vname(3,indxdusthard)='kg m-3                                '
      vname(4,indxdusthard)='  '
      vname(1,indxdustsoft)='dust_soft                             '
      vname(2,indxdustsoft)='Soft component of dust                '
      vname(3,indxdustsoft)='kg m-3                                '
      vname(4,indxdustsoft)='  '

      vname(1,indxpochard)='POC_hard'
      vname(2,indxpochard)='Hard component of POC'
      vname(3,indxpochard)='mmol C/m^3'
      vname(4,indxpochard)='  '
      vname(1,indxpocsoft)='POC_soft'
      vname(2,indxpocsoft)='Soft component of POC'
      vname(3,indxpocsoft)='mmol C/m^3'
      vname(4,indxpocsoft)='  '

      vname(1,indxpcaco3hard)='PCaCO3_hard'
      vname(2,indxpcaco3hard)='Hard component of particulate CaCO3'
      vname(3,indxpcaco3hard)='mmol CaCO3/m^3'
      vname(4,indxpcaco3hard)='  '
      vname(1,indxpcaco3soft)='PCaCO3_soft'
      vname(2,indxpcaco3soft)='Soft component of particulate CaCO3'
      vname(3,indxpcaco3soft)='mmol CaCO3/m^3'
      vname(4,indxpcaco3soft)='  '

      vname(1,indxpsio2hard)='PSiO2_hard'
      vname(2,indxpsio2hard)='Hard component of particulate SiO2'
      vname(3,indxpsio2hard)='mmol Si/m^3'
      vname(4,indxpsio2hard)='  '
      vname(1,indxpsio2soft)='PSiO2_soft'
      vname(2,indxpsio2soft)='Soft component of particulate SiO2'
      vname(3,indxpsio2soft)='mmol Si/m^3'
      vname(4,indxpsio2soft)='  '

      vname(1,indxpironhard)='PFe_hard'
      vname(2,indxpironhard)='Hard component of particulate iron'
      vname(3,indxpironhard)='mmol Fe/m^3'
      vname(4,indxpironhard)='  '
      vname(1,indxpironsoft)='PFe_soft'
      vname(2,indxpironsoft)='Soft component of particulate iron'
      vname(3,indxpironsoft)='mmol Fe/m^3'
      vname(4,indxpironsoft)='  '

# if defined BEC2_DIAG
      vname_bec2_diag_3d(1,indxDUSTHARDREMIN)='dust_hard_remin       '
      vname_bec2_diag_3d(2,indxDUSTHARDREMIN)=
     &       'Remin of Hard component of dust'
      vname_bec2_diag_3d(3,indxDUSTHARDREMIN)='kg m-3 s-1 (check) '
      vname_bec2_diag_3d(4,indxDUSTHARDREMIN)='  '

      vname_bec2_diag_3d(1,indxDUSTSOFTREMIN)='dust_soft_remin        '
      vname_bec2_diag_3d(2,indxDUSTSOFTREMIN)=
     &       'Remin of Soft component of dust'
      vname_bec2_diag_3d(3,indxDUSTSOFTREMIN)='kg m-3 s-1 (check)   '
      vname_bec2_diag_3d(4,indxDUSTSOFTREMIN)='  '

      vname_bec2_diag_3d(1,indxPOCHARDREMIN)='POC_hard_remin'
      vname_bec2_diag_3d(2,indxPOCHARDREMIN)=
     &        'Remin of Hard component of POC'
      vname_bec2_diag_3d(3,indxPOCHARDREMIN)='mmol C m-3 s-1'
      vname_bec2_diag_3d(4,indxPOCHARDREMIN)='  '

      vname_bec2_diag_3d(1,indxPOCSOFTREMIN)='POC_soft_remin'
      vname_bec2_diag_3d(2,indxPOCSOFTREMIN)=
     &        'Remin of Soft component of POC'
      vname_bec2_diag_3d(3,indxPOCSOFTREMIN)='mmol C m-3 s-1'
      vname_bec2_diag_3d(4,indxPOCSOFTREMIN)='  '

      vname_bec2_diag_3d(1,indxCACO3HARDREMIN)='CaCO3_hard_remin'
      vname_bec2_diag_3d(2,indxCACO3HARDREMIN)=
     &       'Remin of Hard component of CaCO3'
      vname_bec2_diag_3d(3,indxCACO3HARDREMIN)='mmol CaCO3 m-3 s-1'
      vname_bec2_diag_3d(4,indxCACO3HARDREMIN)='  '

      vname_bec2_diag_3d(1,indxCACO3SOFTREMIN)='CaCO3_soft_remin'
      vname_bec2_diag_3d(2,indxCACO3SOFTREMIN)=
     &        'Remin of Soft component of CaCO3'
      vname_bec2_diag_3d(3,indxCACO3SOFTREMIN)='mmol CaCO3 m-3 s-1'
      vname_bec2_diag_3d(4,indxCACO3SOFTREMIN)='  '

      vname_bec2_diag_3d(1,indxSIO2HARDREMIN)='SiO2_hard_remin'
      vname_bec2_diag_3d(2,indxSIO2HARDREMIN)=
     &         'Remin of Hard component of SiO2'
      vname_bec2_diag_3d(3,indxSIO2HARDREMIN)='mmol Si m-3 s-1'
      vname_bec2_diag_3d(4,indxSIO2HARDREMIN)='  '

      vname_bec2_diag_3d(1,indxSIO2SOFTREMIN)='SiO2_soft_remin'
      vname_bec2_diag_3d(2,indxSIO2SOFTREMIN)=
     &       'Remin of Soft component of SiO2'
      vname_bec2_diag_3d(3,indxSIO2SOFTREMIN)='mmol Si m-3 s-1'
      vname_bec2_diag_3d(4,indxSIO2SOFTREMIN)='  '

      vname_bec2_diag_3d(1,indxPIRONHARDREMIN)='PFe_hard_remin'
      vname_bec2_diag_3d(2,indxPIRONHARDREMIN)=
     &        'Remin of Hard component of piron'
      vname_bec2_diag_3d(3,indxPIRONHARDREMIN)='mmol Fe m-3 s-1'
      vname_bec2_diag_3d(4,indxPIRONHARDREMIN)='  '

      vname_bec2_diag_3d(1,indxPIRONSOFTREMIN)='PFe_soft_remin'
      vname_bec2_diag_3d(2,indxPIRONSOFTREMIN)=
     &        'Remin of Soft component of piron'
      vname_bec2_diag_3d(3,indxPIRONSOFTREMIN)='mmol Fe m-3 s-1'
      vname_bec2_diag_3d(4,indxPIRONSOFTREMIN)='  '

# endif /* BEC2_DIAG */
# endif /* USE_EXPLICIT_VSINK */


# ifdef BEC2_DIAG

        ! MF: Control BEC2_DIAG Output vars
      wrtavg_bec2_diag_3d=(/.true.,.false.,.false.,.true.,.false.,      !PAR, POC_FLUX_IN, POC_Prod,POC_REMIN, CACO3_FLUX_IN,
     &                     .false.,.false.,.false.,.false.,.false.,       !PCACO3_PROD, CACO3_REMIN,SIO2_FLUX_IN, SIO2_PROD, SIO2_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,       !DUST_FLUX_IN, DUST_REMIN,P_IRON_FLUX_IN, P_IRON_PROD, P_IRON_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,       !GRAZE_SP, GRAZE_DIAT,GRAZE_DIAZ, SP_LOSS, DIAT_LOSS,
     &                     .false.,.false.,.false.,.false.,.false.,       !ZOO_LOSS, SP_AGG, DIAT_AGG,PHOTOC_SP, PHOTOC_DIA,
     &                     .true.,.false.,.true.,.false.,.false.,       !TOT_PROD, DOC_PROD,DOC_REMIN, FE_SCAVENGE, SP_N_LIM,
     &                     .false.,.false.,.false.,.false.,.false.,       !SP_FE_LIM, SP_PO4_LIM,SP_LIGHT_LIM, DIAT_N_LIM, DIAT_FE_LIM,
     &                     .false.,.false.,.false.,.false.,.false.,       !DIAT_PO4_LIM, DIAT_SIO3_LIM,DIAT_LIGHT_LIM, CACO3_PROD, DIAZ_NFIX,
     &                     .false.,.false.,.false.,.false.,.false.,       !DIAZ_LOSS, PHOTOC_DIAZ,DIAZ_PO4_LIM, DIAZ_FE_LIM, DIAZ_LIGHT_LIM,
     &                     .true.,.false.,.false.,.false.,.false.,       !FE_SCAVENGE_RATE, DON_PROD,DON_REMIN, DOFE_PROD, DOFE_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,     !DOP_PROD, DOP_REMIN,DIAT_SI_UPTAKE, PHOTOFE_SP, PHOTOFE_DIA,
     &                     .false.,.false.,.false.,.false.,.false.,    !PHOTOFE_DIAZ, NITRI,DENITRI,SP_NO3_UPTAKE,DIAT_NO3_UPTAKE,
     &                     .false.,.false.,.false.,.false.,.false.,  !DIAZ_NO3_UPTAKE,SP_NH4_UPTAKE, DIAT_NH4_UPTAKE,DIAZ_NH4_UPTAKE,SP_GRAZE_DIC,
     &                     .false.,.false.,.false.,.false.,.false.,  !DIAT_GRAZE_DIC,DIAZ_GRAZE_DIC, SP_LOSS_DIC,DIAT_LOSS_DIC,DIAZ_LOSS_DIC,
     &                     .false.,.false.,.false.,.false.,.false.,  !ZOO_LOSS_DIC, DIAZ_AGG,GRAZE_SP_ZOO, GRAZE_DIAT_ZOO, GRAZE_DIAZ_ZOO,
     &                     .false.,.false.,.false.,.false.,          !SPQCACO3,PHOTOACC_SP,PHOTOACC_DIAT, PHOTOACC_DIAZ,
     &                     .false.,.false.,.false.,.false.,.false.,  !SPC_ZERO, DIATC_ZERO, DIAZC_ZERO, DOC_ZERO, ZOOC_ZERO
     &                     .false.,.false.,.true.,.false.,           !SPCACO3_ZERO, DONR_REMIN, TOTCHL, SP_P_LIM
     &                     .false.,.false.,.false.,.false.,.false.      !DIAT_P_LIM,DIAZ_P_LIM, TOT_PHYTO_C, O2_CONSUMPTION, O2_PRODUCTION
     ! 93 tracer
# ifdef USE_EXPLICIT_VSINK
     ! +10 -> 101 tracer
     &                     ,.false.,.true.,.true.,.false.,.true.,  ! PIRON_HARD_REMIN, CACO3_HARD_REMIN, SIO2_HARD_REMIN, POC_HARD_REMIN,DUST_HARD_REMIN
     &                     .true.,.true.,.true.,.true.,.true.    ! PIRON_SOFT_REMIN, CACO3_SOFT_REMIN, SIO2_SOFT_REMIN, POC_SOFT_REMIN,DUST_SOFT_REMIN
# endif
# ifdef BEC_COCCO
     ! + 18 -> 119 tracer (w/ expl sinking) or 109 tracer (w/ impl sinking)
     &                     ,.true., ! GRAZE_COCCO
     &                     .true.,.true.,.true.,.true.,     ! COCCO_LOSS,COCCO_AGG, PHOTOC_COCCO, COCCO_N_LIM
     &                     .true.,.true.,.true.,.true.,     ! COCCO_PO4_LIM, COCCO_FE_LIM, COCCO_LIGHT_LIM, CACO3_PROD_COCCO
     &                     .true.,.true.,.true.,.true., ! IRON_UPTAKE_COCCO, COCCO_NO3_UPTAKE, COCCO_NH4_UPTAKE, COCCO_GRAZE_DIC
     &                     .true.,.true.,.true.,.true.,    ! COCCO_LOSS_DIC, GRAZE_COCCO_ZOO, QCACO3_COCCO, COCCO_PHOTOACC
     &                     .true.   ! COCCO_P_LIM
# endif
# ifdef Ncycle_SY
     &                     ,.true.,.true.,.true.,.true.,.true.,.true.!AMMOX, NITROX,ANAMMOX,DENITRIF1,DENITRIF2,DENITRIF3
     &                     ,.true.,.true.,.true., .true., .true. !SP_NO2_UPTAKE,DIAT_NO2_UPTAKE, DIAZ_NO2_UPTAKE, N2OAMMOX, N2OSODEN_CONS
     &                     ,.true.,.true.,.true. !N2OAO1_CONS, N2OATM_CONS, N2OSIDEN_CONS
# endif
# ifdef N2O_NEV
     &                     ,.true.,.true. ! n2o_cons_nev, n2o_prod_nev
# endif
     &                     /)

# endif /* BEC2_DIAG */

      tclm_name(iPO4) = 'po4_time'
      tclm_name(iNO3) = 'no3_time'
      tclm_name(iSIO3) = 'sio3_time'
      tclm_name(iO2) = 'o2_time'
      tclm_name(iFE) = 'fe_time'
      tclm_name(iDIC) = 'dic_time'
      tclm_name(iALK) = 'alk_time'
      tclm_name(iNH4) = 'nh4_time'
      tclm_name(iDOC) = 'doc_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOPR) = 'dop_time'
      tclm_name(iDON) = 'don_time'
      tclm_name(iDONR) = 'don_time'
      tclm_name(iDOFE) = 'dofe_time'
      tclm_name(iSPC) = 'spc_time'
      tclm_name(iSPCHL) = 'spchl_time'
      tclm_name(iSPCACO3) = 'spcaco3_time'
      tclm_name(iSPFE) = 'spfe_time'
      tclm_name(iDIATC) = 'diatc_time'
      tclm_name(iDIATCHL) = 'diatchl_time'
      tclm_name(iDIATSI) = 'diatsi_time'
      tclm_name(iDIATFE) = 'diatfe_time'
      tclm_name(iZOOC) = 'zooc_time'
      tclm_name(iDIAZC) = 'diazc_time'
      tclm_name(iDIAZCHL) = 'diazchl_time'
      tclm_name(iDIAZFE) = 'diazfe_time'
# ifdef BEC_COCCO
      tclm_name(iCOCCOC) = 'coccoc_time'
      tclm_name(iCOCCOCHL) = 'coccochl_time'
      tclm_name(iCOCCOCAL) = 'coccocal_time'
      tclm_name(iCOCCOFE) = 'coccofe_time'
      tclm_name(iCAL) = 'cal_time'
# endif
# ifdef Ncycle_SY
      tclm_name(iNO2) = 'no2_time'
      tclm_name(iN2) = 'n2_time'
      tclm_name(iN2O) = 'n2o_time'
# ifdef N2O_TRACER_DECOMP
      tclm_name(iN2O_AO1) = 'n2o_ao1_time'
      tclm_name(iN2O_SIDEN) = 'n2o_siden_time'
      tclm_name(iN2O_SODEN) = 'n2o_soden_time'
      tclm_name(iN2O_ATM) = 'n2o_atm_time'
      tclm_name(iN2_SED) = 'n2_sed_time'
# endif
# endif
# ifdef N2O_NEV
      tclm_name(iN2O_NEV) = 'n2o_nev_time'
# endif
#endif /* SOLVE3D */

      end subroutine init_scalars_bec2
! ----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 - entire module */
      end module bgc
