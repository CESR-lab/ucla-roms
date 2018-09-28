#include "cppdefs.h"
#ifdef STATIONS
c--# define VERBOSE

      subroutine wrt_statn                    ! Write out data into
                                              ! stations NetCDF file
      implicit none
# include "param.h"
# include "ncvars.h"
      real psta(NS)
      integer record, i,j, nst, ierr, start(3),count(3), lvar, lenstr
# ifdef SOLVE3D
     &                                             , itrc, k, tile
      real rsta(N,NS)
# endif
# include "scalars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "mixing.h"
# include "netcdf.inc"
# include "work.h"

#if defined MPI && !defined PARALLEL_FILES
# include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank
      if (mynode>0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

! Create/open averages file; write grid arrays, if so needed, check
! for errors, then set record within the file "record".

      call def_statn (ncidstn, nrecstn, ierr)
      if (ierr /= nf_noerr) goto 99         !!! WARNING: Here it is
      nrecstn=max(nrecstn,1)                  !!! assumed that global
      if (nrpfstn==0) then                  !!! averages record index
        record=nrecstn                        !!! "nrecstn" is already
      else                                    !!! advanced by main.
        record=1+mod(nrecstn-1, nrpfstn)
      endif

      start(1)=1
      count(1)=N
      start(2)=1
      count(2)=nstation
      start(3)=record
      count(3)=1

! Time

      ierr=nf_put_var1_FTYPE (ncidstn, stnTime, record, time)
      if (ierr /= nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) record, vname(1,indxTime)(1:lvar), nf_strerror(ierr)
        goto 99                                           !--> ERROR
      endif


! Free-surface (m).

      if (wrtsta(indxZ)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          psta(nst)=zeta(i,j,knew)
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn,stnZ, start(2),count(2), psta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(*,1) record, vname(1,indxZ)(1:lvar), nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif


! 2D barotropic velocity components in XI- and ETA-directions [m/s]

      if (wrtsta(indxUb)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          psta(nst)=0.5*(ubar(i,j,knew)+ubar(i+1,j,knew))
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn,stnUb, start(2),count(2), psta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          write(*,1) record, vname(1,indxUb)(1:lvar), nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

      if (wrtsta(indxVb)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          psta(nst)=0.5*(vbar(i,j,knew)+vbar(i,j+1,knew))
        enddo
        ierr=nf_put_vara_FTYPE(ncidstn,stnVb, start(2),count(2), psta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          write(*,1) record, vname(1,indxVb)(1:lvar), nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

# ifdef SOLVE3D

! 3D momentum component (m/s) in the XI-direction.

      if (wrtsta(indxU)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=0.5*(u(i,j,k,nstp)+u(i+1,j,k,nstp))
          enddo
        enddo
        ierr=nf_put_vara_FTYPE(ncidstn, stnU, start,count,rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          write(*,1) record, vname(1,indxU)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

! 3D momentum component (m/s) in the ETA-direction.

      if (wrtsta(indxV)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=0.5*(v(i,j,k,nstp)+v(i,j+1,k,nstp))
          enddo
        enddo
        ierr=nf_put_vara_FTYPE(ncidstn, stnV, start,count, rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(*,1) record, vname(1,indxV)(1:lvar), nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

! Vertical velocity (m/s).

      if (wrtsta(indxW)) then
        do tile=0,NSUB_X*NSUB_E-1
          call wvlcty (work,tile)
        enddo
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=work(i,j,k)
          enddo
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn, stnW, start,count, rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          write(*,1) record, vname(1,indxW)(1:lvar), nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

! Tracer variables.

      do itrc=1,NT
        if (wrtsta(itrc+itemp-1)) then
          do nst=1,nstation
            i=ispos(nst)
            j=jspos(nst)
            do k=1,N
              rsta(k,nst)=t(i,j,k,nstp,itrc)
            enddo
          enddo
          ierr=nf_put_vara_FTYPE (ncidstn, stnT(itrc), start,count,
     &                                                       rsta)
          if (ierr/=nf_noerr) then
            lvar=lenstr(vname(1,itrc+itemp-1))
            write(*,1) record, vname(1,itrc+itemp-1)(1:lvar),
     &                                     nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif
      enddo

! Density anomaly.

      if (wrtsta(indxR)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
# ifdef SPLIT_EOS
            rsta(k,nst)=rho1(i,j,k)
# else
            rsta(k,nst)=rho(i,j,k)
# endif
          enddo
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn, stnR, start,count, rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          write(*,1) record,vname(1,indxR)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif
#  ifdef LMD_KPP

!  Depth of planetary boundary layer (m).

      if (wrtsta(indxHbls)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          psta(nst)=hbls(i,j,nstp)
        enddo
        ierr=nf_put_vara_FTYPE(ncidstn,stnHbls, start(2),count(2),psta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxHbls))
          write(*,1) record,vname(1,indxHbls)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif
#  endif /* LMD_KPP */

! Vertical viscosity

      if (wrtsta(indxAkv)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=Akv(i,j,k)
          enddo
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn, stnAkv, start,count, rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(*,1) record,vname(1,indxAkv)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif

! Vertical heat conductivity

      if (wrtsta(indxAkt)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=Akt(i,j,k,itemp)
          enddo
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn, stnAkt, start,count, rsta)
        if (ierr/=nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(*,1) record,vname(1,indxAkt)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif
#  ifdef SALINITY

! Vertical diffusion coefficient for salinity.

      if (wrtsta(indxAks)) then
        do nst=1,nstation
          i=ispos(nst)
          j=jspos(nst)
          do k=1,N
            rsta(k,nst)=Akt(i,j,k,isalt)
          enddo
        enddo
        ierr=nf_put_vara_FTYPE (ncidstn, stnAks, start, count, rsta)
        if (ierr /= nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(*,1) record,vname(1,indxAks)(1:lvar),nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      endif
#  endif /* SALINITY */
# endif /* SOLVE3D */

  1   format(/1x, '### ERROR: wrt_statn :: Cannot write rec =', I5,1x,
     &   'for variable ''', A, ''' into stations NetCDF file:'/23x, A)
      goto 100
  99  if (may_day_flag==0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to the disk to allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close (ncidstn)
      if (nrpfstn>0 .and. record>=nrpfstn) ncidstn=-1
#else
      if (nrpfstn>0 .and. record>=nrpfstn) then
        ierr=nf_close (ncidstn)
        ncidstn=-1
      else
        ierr=nf_sync(ncidstn)
      endif
#endif
      if (ierr/=nf_noerr) then
        write(*,'(/1x,2A/)')     '### ERROR: wrt_statn :: Cannot ',
     &                'synchronize or close stations netCDF file.'
        if (may_day_flag==0) may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                            mynode+1, 2, ocean_grid_comm, ierr)
#endif
      end
#else
      subroutine wrt_statn_empty
      end
#endif /* STATIONS */
