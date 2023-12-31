#include "cppdefs.opt"
#if defined M2NUDGING && !defined M2_FRC_BRY && !defined ANA_SSH

                                      ! Read Sea Surface Height (SSH)
      subroutine get_ssh(ierr_ext)    ! from climatological file.

      use param
      use climat
      use scalars
      use netcdf

      implicit none
      integer :: ierr_ext, ierr, lstr, lenstr, ncdf_read
      real(kind=8) :: dtime
# define SSH_DATA
# include "netcdf.inc"

! Initialization:  Check, whether climatology netCDF file is already
!----------------  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables. Determine whether there is cycling to
! reuse the input data and find cycling period "ssh_cycle", set
! initial cycling index "ssh_ncycle" and record index "ssh_rec".

      ierr=nf90_noerr
      if (iic == ntstart) then

! Note that these three must be initialized here no matter what: not
! doing it causes segmentation fault and uninitialized variable errors
! in do while(ssh_time(itssh),...)  below if something goes wrong with
! reading data, even thought the .and.ierr==nf90_noerr  logic should
! prevent execution of this while. Set both bounding times to
! unrealistacally large negative values to trigger the while-loop.

        itssh=1 ; ssh_time(1)=-1.E+33 ; ssh_time(2)=-1.E+33

        lstr=lenstr(clm_file)
        if (ncidclm == -1) then
          ierr=nf90_open(clm_file(1:lstr), nf_nowrite, ncidclm)
          if (ierr /= nf90_noerr) write(*,'(/1x,4A/)') '### ERROR: ',
     &                    'get_ssh :: can not open netCDF file ''',
     &                                     clm_file(1:lstr), '''.'
        endif
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncidclm, 'ssh_time', ssh_tid)
          if (ierr == nf90_noerr) then
            ierr=nf90_inq_varid(ncidclm, 'zeta', ssh_id)
            if (ierr /= nf90_noerr) then
              ierr=nf90_inq_varid(ncidclm, 'SSH', ssh_id)
            endif
            if (ierr == nf90_noerr) then
              call find_record(ncidclm, ssh_tid,  ssh_rec, ntssh,
     &                         ssh_ncycle, ssh_cycle,  1,1, ierr)

            else
              write(*,'(1x,5A/12x,A/)')     '### ERROR: get_ssh :: ',
     &          'Cannot get netCDF ID for variable "zeta" or "SSH" ',
     &          'from ''', clm_file(1:lstr),'''.', nf90_strerror(ierr)
            endif
          else
            write(*,'(1x,5A/12x,A/)')       '### ERROR: get_ssh :: ',
     &        'Cannot get netCDF ID for variable "ssh_time" from ''',
     &                    clm_file(1:lstr), '''.', nf90_strerror(ierr)
          endif
        endif
      endif

! Read data from the file:    Check if model time is bounded by past
!----- ---- ---- --- -----    and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! data. Repeat until model time falls between the two data times.

      do while(ssh_time(itssh) < time+0.5*dt .and.ierr == nf90_noerr)
        call advance_record(ssh_rec, ntssh,  ssh_ncycle, ssh_cycle,
     &                                                   1,1, ierr)
        if (ierr == nf90_noerr) then
          ierr=nf_get_var1_double(ncidclm, ssh_tid, ssh_rec, dtime)
          if (ierr == nf90_noerr) then
            itssh=min(3-itssh,ntssh)
            ssh_time(itssh)=dtime*day2sec+ssh_cycle*ssh_ncycle
            ierr=ncdf_read(ncidclm, ssh_id, ssh_rec, r_var,
     &                     sshg(START_2D_ARRAY,itssh), 1)
            if (ierr == nf90_noerr) then
              mpi_nonexit_warn write(*,'(7x,A,27x,A,F12.4,2(1x,A,I4))')
     &               'get_ssh :: read SSH', 'ssh_time =', dtime,
     &               'rec =', ssh_rec MYID
              if (ntssh == 1) return
            else
              write(*,2)    'SSH', ssh_rec,  nf90_strerror(ierr)
            endif
          else
            write(*,2) 'ssh_time', ssh_rec,  nf90_strerror(ierr)
          endif
        else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &   '### ERROR: get_ssh :: Requested time record ',  ssh_rec,
     &   'exeeds the last record',  ntssh, 'available in netCDF ',
     &   'file ''',  clm_file(1:lstr),  '''',  'tdays = ',  tdays,
     &   '  but the last available ssh_time =',
     &                                    ssh_time(itssh)*sec2day
        endif
      enddo
      if (ierr_ext == nf90_noerr) ierr_ext=ierr

  2   format(/1x,'### ERROR: get_ssh :: Cannot read variable ''',
     &                              A, ''' for rec =', I4/12x,A)
      end



      subroutine set_ssh_tile(istr,iend,jstr,jend, ierr)

! Set-up sea surface height climatology for current tile.

      use param
      use climat
      use scalars

      implicit none
      integer istr,iend,jstr,jend, ierr, i,j, it1,it2
      real tmid, cff1,cff2, cff
# define SSH_DATA

# include "compute_extended_bounds.h"

! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.

      if (ntssh > 1) then                  !--> interpolate ssh in time
        tmid=time+0.5*dt
        it2=itssh ; cff1=ssh_time(it2)-tmid
        it1=3-it2 ; cff2=tmid-ssh_time(it1)
        if (cff1 < dt) synchro_flag=.true.
        cff=1./(cff1+cff2)
        cff1=cff1*cff ; cff2=cff2*cff

        do j=jstrR,jendR
          do i=istrR,iendR
            ssh(i,j)=cff1*sshg(i,j,it1)+cff2*sshg(i,j,it2)
          enddo
        enddo

        if (cff1 < 0. .or. cff2 < 0.) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')      '### WARNING: ',
     &       'set_ssh_tile :: Current model time is outside bounds ',
     &     'of ''ssh_time''.', 'ssh_tstart =', ssh_time(it2)*sec2day,
     &                         'tdays =',      tdays,
     &                         'ssh_tend =',   ssh_time(it2)*sec2day
          endif
          ierr=ierr+1
        endif
      elseif (iic == ntstart) then        ! load time-invariant
        do j=jstrR,jendR                  ! sea-surface height
          do i=istrR,iendR
            ssh(i,j)=sshg(i,j,itssh)
          enddo
        enddo
      endif
      end
#else
      subroutine get_ssh_empty
      end
#endif
