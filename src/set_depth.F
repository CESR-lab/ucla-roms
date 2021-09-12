#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine set_depth(tile)

      use param

      implicit none
      integer tile
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_depth_tile(istr,iend,jstr,jend)
      end

      subroutine set_depth_tile(istr,iend,jstr,jend)

! Compute evolving z-coordinates "z_r" and "z_w" and the heights
! of vertical grid boxes "Hz" from fast-time-averaged free surface,
! topography h(i,j), and vertical coordinate transformation formula.

      use param
# ifdef WEC
      use wec_frc
# endif
      use coupling

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff_r,cff1_r, cff_w,cff1_w, ds

# include "scoord.h"
# include "scalars.h"

# ifndef VERT_COORD_TYPE_SM09
      real cff2_r,cff2_w, z_r0,z_w0
# endif
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# include "compute_auxiliary_bounds.h"

                                       ! During initialization or
      if (iic == 0) then               ! restart compute barotropic
        do j=jstrR,jendR               ! fluxes and place them into
          do i=istrR,iendR             ! arrays for fast-time-averaged
            hinv(i,j)=1./(h(i,j)+hc)   ! fluxes. Later during the first
          enddo                        ! procedure uses these DU_avg1,
        enddo                          ! DV_avg1 as input to compute
                                       ! ubar,vbar as it does so during
        do j=jstrR,jendR               ! all subsequent steps.
          do i=istr,iendR
            DU_avg1(i,j)=0.5*( h(i-1,j)+h(i,j)+zeta(i,j,1)
     &            +zeta(i-1,j,1) )*dn_u(i,j)*( ubar(i,j,1)
# ifdef WEC
     &                                        + ust2d(i,j)
# endif
     &                                                   )
          enddo
        enddo

        do j=jstr,jendR
          do i=istrR,iendR
            DV_avg1(i,j)=0.5*( h(i,j)+h(i,j-1)+zeta(i,j,1)
     &            +zeta(i,j-1,1) )*dm_v(i,j)*( vbar(i,j,1)
# ifdef WEC
     &                                        + vst2d(i,j)
# endif
     &                                                   )
          enddo
        enddo

# ifdef OBC_CHECK
        print *, 'set_depth_tile:'
        print *, 'DU_avg1(istr,10)=',DU_avg1(istr,10)
# endif

      endif

      ds=1.D0/dble(N)
      do j=jstrR,jendR
        do i=istrR,iendR
          z_w(i,j,0)=-h(i,j)
        enddo

        do k=1,N,+1   !--> irreversible because of recursion in Hz

          cff_w=hc*ds* dble(k-N)
          cff_r=hc*ds*(dble(k-N)-0.5)

          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)

          do i=istrR,iendR

            z_w(i,j,k)=zeta(i,j,knew) +(zeta(i,j,knew)+h(i,j))
     &                        *(cff_w+cff1_w*h(i,j))*hinv(i,j)

            z_r(i,j,k)=zeta(i,j,knew) +(zeta(i,j,knew)+h(i,j))
     &                        *(cff_r+cff1_r*h(i,j))*hinv(i,j)

            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
          enddo
        enddo
      enddo
# ifdef EXCHANGE
      if (iic == 0) then
        call exchange2d_tile(istr,iend,jstr,jend, hinv)
      endif
      call exchange_3_tile(istr,iend,jstr,jend, z_w,N+1, z_r,N, Hz,N)
# endif

#if defined NHMG || defined NONTRAD_COR
!     Compute slopes
      do k = 1, N
         do j = 0,Mm+1
            do i = 0,Lm+1
               dzdxi (i,j,k)=0.5*(z_r(i+1,j,k)-z_r(i-1,j,k))*pn(i,j)
               dzdeta(i,j,k)=0.5*(z_r(i,j+1,k)-z_r(i,j-1,k))*pm(i,j)
            enddo
         enddo
      enddo
      ! z_r is buffer filled in the interior halos but not in the outer halos
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
         do k = 1, N
            do j = 0,Mm+1
               dzdxi (0,j,k) = 0.
               dzdxi (1,j,k) = 0.
               dzdeta(0,j,k) = 0.
               dzdeta(1,j,k) = 0.
            enddo
         enddo
       endif
      if (EASTERN_EDGE) then
         do k = 1, N
            do j = 0,Mm+1
               dzdxi (Lm  ,j,k) = 0.
               dzdxi (Lm+1,j,k) = 0.
               dzdeta(Lm  ,j,k) = 0.
               dzdeta(Lm+1,j,k) = 0.
            enddo
         enddo
       endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
         do k = 1, N
            do i = 0,Lm+1
               dzdxi (i,0,k) = 0.
               dzdxi (i,1,k) = 0.
               dzdeta(i,0,k) = 0.
               dzdeta(i,1,k) = 0.
            enddo
         enddo
       endif
      if (NORTHERN_EDGE) then
         do k = 1, N
            do i = 0,Lm+1
               dzdxi (i,Mm  ,k) = 0.
               dzdxi (i,Mm+1,k) = 0.
               dzdeta(i,Mm  ,k) = 0.
               dzdeta(i,Mm+1,k) = 0.
            enddo
         enddo
       endif
# endif

#endif /* NHMG */


      end



      subroutine set_HUV(tile)

      use param

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call set_HUV_tile (istr,iend,jstr,jend)
      end

      subroutine set_HUV_tile (istr,iend,jstr,jend)

      use param
# ifdef WEC
      use wec_frc
# endif

      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            FlxU(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)

     &                                       *( u(i,j,k,nrhs)
# ifdef WEC 
     &                                            +ust(i,j,k)
# endif
     &                                                      )
! JM make sure that ust,vst are properly masked in wec_frc and remove
# if defined WEC && defined MASKING
            FlxU(i,j,k)=FlxU(i,j,k)*umask(i,j)
# endif
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
            FlxV(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
! Daniel add/edit for WEC
!     &                                         *v(i,j,k,nrhs)
     &                                       *( v(i,j,k,nrhs)
# ifdef WEC 
     &                                           + vst(i,j,k)
# endif
     &                                                      )
# if defined WEC && defined MASKING
            FlxV(i,j,k)=FlxV(i,j,k)*vmask(i,j)
# endif
          enddo
        enddo
      enddo
# ifdef EXCHANGE
      call exchange_2_tile (istr,iend,jstr,jend, FlxU,N, FlxV,N)
# endif
      end




      subroutine set_HUV1(tile)

      use param

      implicit none
      integer tile
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_HUV1_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2))
      end

      subroutine set_HUV1_tile (istr,iend,jstr,jend, DC,FC)

      use param
# ifdef WEC
      use wec_frc
# endif
      use coupling

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: DC,FC
# include "grid.h"
# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      real cff
      integer kbak
#  include "ocean2d.h"
# endif
# include "ocean3d.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      kbak=knew-1
      if (kbak < 1) kbak=4
# endif

      do j=jstrR,jendR ! big j-loop

        do i=istr,iendR
          DC(i,N)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))*dn_u(i,j)
          DC(i,0)=DC(i,N)
          FC(i,0)=DC(i,N)*u(i,j,N,nnew)
        enddo
        do k=N-1,1,-1
          do i=istr,iendR
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo

c-->#undef EXTRAP_BAR_FLUXES  /* <-- for testing purposes only */

# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
        if (FIRST_TIME_STEP) then
#  endif
          do i=istr,iendR
            FC(i,0)=(FC(i,0)-DU_avg1(i,j))/DC(i,0)
#  ifdef WEC 
     &                                 +ust2d(i,j)
#  endif
          enddo

#  ifdef EXTRAP_BAR_FLUXES
! PAC23 verified setting with theoretical alpha_max=1.0877.
! This should match the simlar-labelled setting in step2d_FB.F.
! All the other settings below  are to be considered as
! experimental.

#   define NOW 3.63D0
#   define MID 4.47D0
#   define BAK 2.05D0



/*
#   define CURV 2.80D0
#   define NOW (1.515D0+CURV)
#   define MID (2.D0*CURV)
#   define BAK (CURV-0.515D0)
*/

/*
#   define NOW 4.D0
#   define MID 5.D0
#   define BAK 2.D0
*/

/*
#   define NOW 3.75D0
#   define MID 4.375D0
#   define BAK 1.5D0
*/

/*
#   define NOW 2.D0
#   define MID 1.D0
#   define BAK 0.D0
*/

!  DC(i,0) is the barotropic surface H*dy = dy*integral(Hz)
!  FC(i,0) is the barotropic flux dy*integral(Hz*u)

        else ! not first time step
          do i=istr,iendR
#   ifdef KEEP_CORIOLIS
            cff=1./DC(i,0)
            ubar(i,j,kbak)=cff*( NOW*DU_avg1(i,j)
     &                          -MID*DU_avg2(i,j)
     &                       +BAK*DU_avg_bak(i,j)
     &                                          )
            FC(i,0)=cff*FC(i,0)-ubar(i,j,kbak)
#   else
            FC(i,0)=( FC(i,0) -NOW*DU_avg1(i,j)
     &                        +MID*DU_avg2(i,j)
     &                     -BAK*DU_avg_bak(i,j)
     &               )/DC(i,0)
#   endif
          enddo
        endif
#  endif
  !! now FC(i,0) = ubar
# else 
        do i=istr,iendR
          FC(i,0)=(FC(i,0)-DU_avg2(i,j))/DC(i,0) !! 2 corrector of ubar? JM
#  ifdef WEC
     &                               +ust2d(i,j)
#  endif
        enddo
# endif

# ifdef OBC_CHECK
        if(j==10) then
          print *, 'set_HUV1: start'
          print *, 'nnew=', nnew
          print *, '   u(istr  ,10,10,nnew)=', u(istr,10,10,nnew)
          print *, '   FC(istr,0)=', FC(istr,0)
          print *, '   u(iend+1,10,10,nnew)=', u(iend+1,10,10,nnew)
          print *, '   FC(iend+1,0)=', FC(iend+1,0)
          print *, '   ubar(istr,10,kbak)=',ubar(istr,10,kbak)
        endif ! debug
# endif

        do k=1,N
          do i=istr,iendR
!!JM make sure that ust, vst is masked in wec_frc and fix this!!
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0)) !! removing ubar from u
# ifdef MASKING
     &                                *umask(i,j)
#  ifdef WEC
     &               +ust(i,j,k)*(umask(i,j)-1.0)
#  endif
# endif
            FlxU(i,j,k)=DC(i,k)*( u(i,j,k,nnew)
# ifdef WEC
     &                                +ust(i,j,k)
# endif
     &                                          )           
          enddo
        enddo

# ifdef OBC_CHECK
        if(j==10) then
          print *, 'set_HUV1: end'
          print *, 'nnew=', nnew
          print *, 'FlxU(istr  ,10,10)=', FlxU(istr,10,10)
          print *, '   u(istr  ,10,10,nnew)=', u(istr,10,10,nnew)
          print *, 'FlxU(iend+1,10,10)=', FlxU(iend+1,10,10)
          print *, '   u(iend+1,10,10,nnew)=', u(iend+1,10,10,nnew)
        endif ! debug
# endif
        if (j >= jstr) then
          do i=istrR,iendR
            DC(i,N)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))*dm_v(i,j)
            DC(i,0)=DC(i,N)
            FC(i,0)=DC(i,N)*v(i,j,N,nnew)
          enddo
          do k=N-1,1,-1
            do i=istrR,iendR
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
          if (FIRST_TIME_STEP) then
#  endif
            do i=istrR,iendR
              FC(i,0)=(FC(i,0)-DV_avg1(i,j))/DC(i,0)
#  ifdef WEC 
     &                                   +vst2d(i,j)
#  endif
            enddo

#  ifdef EXTRAP_BAR_FLUXES
          else
            do i=istrR,iendR
#   ifdef KEEP_CORIOLIS
              cff=1./DC(i,0)
              vbar(i,j,kbak)=cff*( NOW*DV_avg1(i,j)
     &                            -MID*DV_avg2(i,j)
     &                         +BAK*DV_avg_bak(i,j)
     &                                            )
              FC(i,0)=cff*FC(i,0)-vbar(i,j,kbak)
#   else
              FC(i,0)=( FC(i,0) -NOW*DV_avg1(i,j)
     &                          +MID*DV_avg2(i,j)
     &                       -BAK*DV_avg_bak(i,j)
     &                 )/DC(i,0)
! Daniel add --> not sure if WEC should be in here. DevinD also not sure but suspect so.
#ifdef WEC
     &                   +vst2d(i,j)
#endif
#   endif
            enddo
          endif
#  endif
# else
          do i=istrR,iendR
            FC(i,0)=(FC(i,0)-DV_avg2(i,j))/DC(i,0)
#  ifdef WEC 
     &                                 +vst2d(i,j)
#  endif
          enddo
# endif
          do k=1,N
            do i=istrR,iendR
!! JM same thing as above for ust. (masking in wec_frc)
              v(i,j,k,nnew)=(v(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                  *vmask(i,j)
#  ifdef WEC
     &                 +vst(i,j,k)*(vmask(i,j)-1.0)
#  endif
# endif
              FlxV(i,j,k)=DC(i,k)*( v(i,j,k,nnew)
# ifdef WEC 
     &                                +vst(i,j,k)
# endif
     &                                          )
            enddo
          enddo
        endif
      enddo  ! <-- j-loop

# ifdef OBC_CHECK
      print *, 'FlxV(10,jstr,10)=', FlxV(10,jstr,10)
      print *, '   v(10,jstr,10,nnew)=', v(10,jstr,10,nnew)
      print *, 'FlxV(10,jend+1,10)=', FlxV(10,jend+1,10)
      print *, '   v(10,jend+1,10,nnew)=', v(10,jend+1,10,nnew)
# endif

# ifdef EXCHANGE
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      call exchange_3_tile(istr,iend,jstr,jend,     FlxU, N,
     &                                u(START_2D_ARRAY,1,nnew), N,
     &                               ubar(START_2D_ARRAY,kbak), 1)

      call exchange_3_tile(istr,iend,jstr,jend,     FlxV, N,
     &                                v(START_2D_ARRAY,1,nnew), N,
     &                               vbar(START_2D_ARRAY,kbak), 1)
#  else
      call exchange_4_tile(istr,iend,jstr,jend, FlxU, N,  FlxV, N,
     &                                u(START_2D_ARRAY,1,nnew), N,
     &                                v(START_2D_ARRAY,1,nnew), N)
#  endif
# endif
      end


      subroutine check_set_HUV1_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param

      implicit none
      integer ierr, is,ie, lenstr
# include "strings.h"
      integer lstr
      character(len=16) tmpstr
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+12
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<set_depth.F>'
# ifdef NOW
      write(tmpstr,'(F16.12)') NOW ; lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+3+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='NOW='/ /tmpstr(1:lstr)
# endif
# ifdef MID
      write(tmpstr,'(F16.12)') MID ; lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+3+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='MID='/ /tmpstr(1:lstr)
# endif
# ifdef BAK
      write(tmpstr,'(F16.12)') BAK ; lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr) == '0' .and. lstr > 3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+3+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='BAK='/ /tmpstr(1:lstr)

! Note that the value of (NOW-MID+BAK-1.D0)/BAK should match
! the coefficient in front of DU_avg2 when computing DU_avg_bak
! in step2d_FB.F. Put it into the signature for cross-checking.

      write(tmpstr,'(F16.12)') (NOW-MID+BAK-1.D0)/BAK
      lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr) == '0' .and. lstr > 3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+9+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='(N-M+B)/B='/ /tmpstr(1:lstr)
# endif
      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine set_depth_empty
      end
#endif /* SOLVE3D */
