!
! $Id: init_arrays_floats.F,v 2.1 2004/04/09 21:23:33 capet Exp $
!
#include "cppdefs.h"
#ifdef FLOATS
      subroutine init_arrays_floats

      implicit none
#include "param.h"
#include "floats.h" 
#include "scalars.h"

      integer iflt, j, i
!
!
!  Initialize floats trajectory array.
!
      do iflt=1,Mfloats
        fltgrd(iflt)=-1 ! -1 is the dead float value. 
        nrecvel(iflt)=0

        do j=igrd,izgrd  
           Tinfo(j,iflt)=spval
        enddo
        
        do j=1,NDIAGS
        trackaux(j,iflt)=0.
        enddo

        do j=0,NFT
          do i=1,NFV
            track(i,j,iflt)=spval
          enddo
        enddo

      enddo

# ifdef AGRIF
      do j=1,maxgrids   ! consistent with initialization of nfp1
        floattindex(j)=0
      enddo
# endif
      return
      end

#else
      subroutine init_arrays_floats_empty
      return
      end

#endif /* FLOATS */

