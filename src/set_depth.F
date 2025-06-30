#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine set_depth(tile)

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call set_depth_tile(istr,iend,jstr,jend)
      end

      subroutine set_depth_tile(istr,iend,jstr,jend)

! Compute evolving z-coordinates "z_r" and "z_w" and the heights
! of vertical grid boxes "Hz" from fast-time-averaged free surface,
! topography h(i,j), and vertical coordinate transformation formula.

      use param
      use wec_frc
      use coupling
      use grid
      use ocean_vars
      use scalars
      use scoord

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff_r,cff1_r, cff_w,cff1_w, ds

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

      endif  ! first time

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
        call exchange_xxx(hinv)
      endif
      call exchange_xxx(z_w,z_r,Hz)
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



      subroutine set_HUV

      use param
      use hidden_mpi_vars

      implicit none
      integer,save:: tile=0
# include "compute_tile_bounds.h"
      call set_HUV_tile (istr,iend,jstr,jend)
      end

      subroutine set_HUV_tile (istr,iend,jstr,jend)

      use param
      use wec_frc
      use grid
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k

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
      call exchange_xxx(FlxU,FlxV)
# endif
      end




      subroutine set_HUV1(tile)

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call set_HUV1_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2))
      end

      subroutine set_HUV1_tile (istr,iend,jstr,jend, DC,FC)

      use param
      use wec_frc
      use coupling
      use grid
      use ocean_vars
      use scalars
      use dimensions
      use diagnostics

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: DC,FC

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      real cff
      integer kbak
# endif

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
          !! JM: Very first time step only

#  ifdef EXTRAP_BAR_FLUXES
! PAC23 verified setting with theoretical alpha_max=1.0877.
! This should match the simlar-labelled setting in step2d_FB.F.
! All the other settings below  are to be considered as
! experimental.

#   define NOW 3.63D0
#   define MID 4.47D0
#   define BAK 2.05D0


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
            FC(i,0)=cff*FC(i,0)-ubar(i,j,kbak)  ! now FC(i,0) is mismatch in[m/s]
#   else
            FC(i,0)=( FC(i,0) -NOW*DU_avg1(i,j)
     &                        +MID*DU_avg2(i,j)
     &                     -BAK*DU_avg_bak(i,j)
     &               )/DC(i,0)
#   endif
          enddo
        endif
#  endif
     !! now FC(i,0) => mismatch between  int(U*dz)/int(dz) and ubar
# else
        do i=istr,iendR
          FC(i,0)=(FC(i,0)-DU_avg2(i,j))/DC(i,0) ! FC(i,0) is mismatch in [m/s]
#  ifdef WEC
     &                               +ust2d(i,j)
#  endif
        enddo
# endif

        do k=1,N
          do i=istr,iendR
!!JM make sure that ust, vst are masked in wec_frc and fix this!!
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0)) !! removing ubar mismatch from u
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
          do k=1,nz
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

# ifdef EXCHANGE
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      call exchange_xxx(FlxU,u(:,:,:,nnew),ubar)
      call exchange_xxx(FlxV,v(:,:,:,nnew),vbar)
#  else
      call exchange_xxx(FlxU,FlxV,u(:,:,:,nnew),v(:,:,:,nnew))
#  endif
# endif
      end


      subroutine check_set_HUV1_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param
      use strings

      implicit none
      integer ierr, is,ie, lenstr

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
      cpps(is:ie)='(N-M+B-1)/B='/ /tmpstr(1:lstr)
# endif
      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine set_depth_empty
      end
#endif /* SOLVE3D */
