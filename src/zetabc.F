#include "cppdefs.opt"

      subroutine zetabc_tile (istr,iend,jstr,jend, zeta_new)

! Set boundary conditions for free surface zeta_new

      use param
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real zeta_new(PRIVATE_2D_SCRATCH_ARRAY), cx, cff,dft,dfx
      real, parameter :: eps=1.E-33

#include "compute_auxiliary_bounds.h"


c--#undef OBC_M2FLATHER
c--#define OBC_Z2ORLANSKI


#ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
# if defined OBC_WEST && defined OBC_M2FLATHER
        do j=jstrV-1,jend
          cx=dtfast*pm(istr,j)*sqrt(g*h(istr,j))
          zeta_new(istr-1,j)=( (1.-cx)*zeta(istr-1,j,kstp)
     &                               +cx*zeta(istr,j,kstp) )
#  ifdef MASKING
     &                                      *rmask(istr-1,j)
#  endif
        enddo
# elif defined OBC_WEST && defined OBC_Z2ORLANSKI
        do j=jstrV-1,jend                             ! Western edge
          dft=zeta(istr,j,kstp)-zeta_new(istr,j)      ! Radiation BC
          dfx=zeta_new(istr,j)-zeta_new(istr+1,j)

          cx=dft*dfx
          if (cx < 0.) cx=0.

          cff=max(dfx*dfx, eps)
          zeta_new(istr-1,j)=( cff*zeta(istr-1,j,kstp)
     &                  +cx*zeta_new(istr,j) )/(cff+cx)

#  ifdef MASKING
     &                                 *rmask(istr-1,j)
#  endif
        enddo
# else
        do j=jstrV-1,jend                              ! Western edge
          zeta_new(istr-1,j)=zeta_new(istr,j)          ! zero gradient BC
#  ifdef MASKING
     &                        *rmask(istr-1,j)
#  endif
        enddo
# endif
      endif

c*******
c*******
# ifdef ISWAKE_XXX
c--#  define OBC_Z2ORLANSKI
# endif
c*******
c*******

      if (EASTERN_EDGE) then
# if defined OBC_EAST && defined OBC_M2FLATHER
        do j=jstrV-1,jend
          cx=dtfast*pm(iend,j)*sqrt(g*h(iend,j))
          zeta_new(iend+1,j)=( (1.-cx)*zeta(iend+1,j,kstp)
     &                               +cx*zeta(iend,j,kstp) )

#  ifdef MASKING
     &                                      *rmask(iend+1,j)
#  endif
        enddo

# elif defined OBC_EAST && defined OBC_Z2ORLANSKI
        do j=jstrV-1,jend                             ! Eastern edge
          dft=zeta(iend,j,kstp)-zeta_new(iend,j)      ! Radiation BC
          dfx=zeta(iend,j,kstp)-zeta(iend-1,j,kstp)
          dfxi=zeta_new(iend,j)-zeta_new(iend-1,j)

          cx=dft*dfx ; cxi=dft*dfxi ;  cff=dfx*dfx
          if (cx < 0.) then
            cx=0. ; cxi=0.
          elseif (cx > c0*cff) then
            cxi=

            write(*,*) 'supercourant'

            c

XXXX TO BE DEVELOPED


          endif

          zeta_new(iend+1,j)=( cffi*zeta(iend+1,j,kstp)-cx(zeta(iend+1,j,kstp)
     &                  +cxi*zeta_new(iend,j) )/(cffi+cxi)



c          cx=dft*dfx
c          if (cx < 0.) cx=0.

          cff=max(dfx*dfx, eps)
          zeta_new(iend+1,j)=( cff*zeta(iend+1,j,kstp)
     &                  +cx*zeta_new(iend,j) )/(cff+cx)
#  ifdef MASKING
     &                                 *rmask(iend+1,j)
#  endif
        enddo
# else
        do j=jstrV-1,jend                             ! Eastern edge
          zeta_new(iend+1,j)=zeta_new(iend,j)         ! gradient BC
#  ifdef MASKING
     &                        *rmask(iend+1,j)
#  endif
        enddo
# endif
      endif
#endif        /* !EW_PERIODIC */


#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
# if defined OBC_SOUTH && defined OBC_M2FLATHER
        do i=istrU-1,iend
          cx=dtfast*pn(i,jstr)*sqrt(g*h(i,jstr))     ! Southern edge
          zeta_new(i,jstr-1)=( (1.-cx)*zeta(i,jstr-1,kstp)
     &                               +cx*zeta(i,jstr,kstp) )
#  ifdef MASKING
     &                                      *rmask(i,jstr-1)
#  endif
        enddo

# elif defined OBC_SOUTH && defined OBC_Z2ORLANSKI
        do i=istrU-1,iend                             ! Southern edge
          dft=zeta(i,jstr,kstp)-zeta_new(i,jstr)      ! Radiation BC
          dfx=zeta_new(i,jstr)-zeta_new(i,jstr+1)

          cx=dft*dfx
          if (cx < 0.) cx=0.

          cff=max(dfx*dfx, eps)

          zeta_new(i,jstr-1)=( cff*zeta(i,jstr-1,kstp)
     &                  +cx*zeta_new(i,jstr) )/(cff+cx)

#  ifdef MASKING
     &                                 *rmask(i,jstr-1)
#  endif
        enddo
# else
        do i=istrU-1,iend                            ! Southern edge
          zeta_new(i,jstr-1)=zeta_new(i,jstr)        ! gradient BC
#  ifdef MASKING
     &                       *rmask(i,jstr-1)
#  endif
        enddo
# endif  /* OBC_SOUTH */
      endif


      if (NORTHERN_EDGE) then
# if defined OBC_NORTH && defined OBC_M2FLATHER
        do i=istrU-1,iend
          cx=dtfast*pn(i,jend)*sqrt(g*h(i,jend))    ! Northern edge
          zeta_new(i,jend+1)=( (1.-cx)*zeta(i,jend+1,kstp)
     &                               +cx*zeta(i,jend,kstp) )
#  ifdef MASKING
     &                                      *rmask(i,jend+1)
#  endif
        enddo

# elif defined OBC_EAST && defined OBC_Z2ORLANSKI
        do i=istrU-1,iend                             ! Northern edge
          dft=zeta(i,jend,kstp)-zeta_new(i,jend)      ! Radiation BC
          dfx=zeta_new(i,jend)-zeta_new(i,jend-1)

          cx=dft*dfx
          if (cx < 0.) cx=0.

          cff=max(dfx*dfx, eps)
          zeta_new(i,jend+1)=( cff*zeta(i,jend+1,kstp)
     &                  +cx*zeta_new(i,jend) )/(cff+cx)
#  ifdef MASKING
     &                                 *rmask(i,jend+1)
#  endif
        enddo

# else
        do i=istrU-1,iend                           ! Northern edge
          zeta_new(i,jend+1)=zeta_new(i,jend)       ! gradient BC
#  ifdef MASKING
     &                       *rmask(i,jend+1)
#  endif
        enddo
# endif   /* OBC_NORTH */
      endif
#endif       /* !NS_PERIODIC */

! Corners between adjacent open boundaries

#ifndef EW_PERIODIC
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
        zeta_new(istr-1,jstr-1)=0.5*( zeta_new(istr  ,jstr-1)
     &                               +zeta_new(istr-1,jstr  ))
      endif
      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
        zeta_new(iend+1,jstr-1)=0.5*( zeta_new(iend  ,jstr-1)
     &                               +zeta_new(iend+1,jstr  ))
      endif
      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
        zeta_new(istr-1,jend+1)=0.5*( zeta_new(istr  ,jend+1)
     &                               +zeta_new(istr-1,jend  ))
      endif
      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
        zeta_new(iend+1,jend+1)=0.5*( zeta_new(iend  ,jend+1)
     &                               +zeta_new(iend+1,jend  ))
      endif
# endif
#endif
      end

