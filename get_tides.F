#include "cppdefs.h"
#if defined SSH_TIDES || defined UV_TIDES
                                     ! Read tidal forcing parameters
      subroutine get_tides (ierr)    ! from netCDF file. Also convert
      implicit none                  ! their units where needed.

# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "tides.h"
# include "ncvars.h"
# include "netcdf.inc"
      integer ierr, ierr_count,  ifrc, ncid, varid,
     &        itide, i,j,  lfrc, lenstr,  ncdf_read
      character*(max_name_size) frcname

      if (iic==ntstart) then
        ifrc=1                 ! since tidal forcing data is read
 100    ncid=ncfrc(ifrc)     ! only once and does not to be re-read
        frcname=frcfile(ifrc)  ! later, there is no need to memorize
        lfrc=lenstr(frcname)   ! netcdf ID for tidal file (unlike for
        ierr=nf_noerr          ! all other forces).
        ierr_count=0

        if (ncid==-1) then
          ierr=nf_open (frcname(1:lfrc), nf_nowrite, ncid)
          if (ierr == nf_noerr) then
            ncfrc(ifrc)=ncid
          else
            write(*,'(/1x,4A/)') '### ERROR in get_tides: Can not ',
     &                'open netCDF file ''', frcname(1:lfrc), '''.',
     &                                            nf_strerror(ierr)
          endif
        endif
        if (ierr == nf_noerr) then
          do itide=1,Ntides
            ierr=nf_inq_varid (ncid, 'tide_period', varid)
            if (ierr == nf_noerr) then
              ierr=nf_get_var1_FTYPE (ncid, varid, itide,
     &                                           Tperiod(itide))
              if (ierr == nf_noerr) then
                write(*,'(5x,2A,I3,1x,3A)')     'get_tides :: read ',
     &                  '''tide_period'' for mode', itide, 'from ''',
     &                                        frcname(1:lfrc), '''.'
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot read ',
     &           'variable ''tide_period'' from ''', frcname(1:lfrc),
     &                                     '''.',  nf_strerror(ierr)
              endif
# ifdef SSH_TIDES

! Read tidal elevation amplitude and phase.

              ierr=nf_inq_varid (ncid, 'tide_Eamp', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid, varid, itide, r_var,
     &                         SSH_Tamp(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(1x,2A,I3,1x,3A)')   'get_tides :: read ',
     &                    '''tide_Eamp'' for mode', itide, 'from ''',
     &                                         frcname(1:lfrc),'''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Eamp'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &                'netCDF ID  variable ''tide_Eamp'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif

              ierr=nf_inq_varid (ncid, 'tide_Ephase', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid, varid, itide, r_var,
     &                         SSH_Tphase(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(1x,2A,I3,1x,3A)')   'get_tides :: read ',
     &                  '''tide_Ephase'' for mode', itide, 'from ''',
     &                                        frcname(1:lfrc), '''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Ephase'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &              'netCDF ID  variable ''tide_Ephase'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif
# endif
# ifdef UV_TIDES

! Read tidal currents angle, phase, major and minor ellipse axis.

              ierr=nf_inq_varid (ncid, 'tide_Cangle', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid, varid, itide, r_var,
     &                         UV_Tangle(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(5x,2A,I3,1x,3A)') 'get_tides :: read ',
     &                '''tide_Cangle'' for mode', itide, 'from ''',
     &                                      frcname(1:lfrc), '''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Cangle'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &              'netCDF ID  variable ''tide_Cangle'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif

              ierr=nf_inq_varid (ncid, 'tide_Cphase', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid, varid, itide, r_var,
     &                         UV_Tphase(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(1x,2A,I3,1x,3A)')  'get_tides :: read ',
     &                 '''tide_Cphase'' for mode', itide, 'from ''',
     &                                       frcname(1:lfrc), '''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Cphase'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &              'netCDF ID  variable ''tide_Cphase'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif

              ierr=nf_inq_varid (ncid, 'tide_Cmax', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid, varid, itide, r_var,
     &                         UV_Tmajor(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(5x,2A,I3,1x,3A)')   'get_tides :: read ',
     &                    '''tide_Cmax'' for mode', itide, 'from ''',
     &                                        frcname(1:lfrc), '''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Cmax'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &                'netCDF ID  variable ''tide_Cmax'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif

              ierr=nf_inq_varid (ncid, 'tide_Cmin', varid)
              if (ierr == nf_noerr) then
                ierr=ncdf_read(ncid,  varid,  itide,  r_var,
     &                         UV_Tminor(START_2D_ARRAY,itide), 1)
                if (ierr == nf_noerr) then
                  write(*,'(5x,2A,I3,1x,3A)')   'get_tides :: read ',
     &                    '''tide_Cmin'' for mode', itide, 'from ''',
     &                                        frcname(1:lfrc), '''.'
                else
                  write(*,'(/1x,4A/8x,A/)')     '### ERROR: Cannot ',
     &                    'read variable ''tide_Cmin'' from ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                  ierr_count=ierr_count+1
                endif
              else
                write(*,'(/1x,4A/8x,A/)')  '### ERROR: Cannot find ',
     &              'netCDF ID  variable ''tide_Cmin'' in file ''',
     &                    frcname(1:lfrc), '''.',  nf_strerror(ierr)
                ierr_count=ierr_count+1
              endif
# endif

! Convert tidal phase and angle from degrees to radians.

              Tperiod(itide)=Tperiod(itide)*3600.
              do j=0,Mm+1
                do i=0,Lm+1
# ifdef SSH_TIDES
                SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)*deg2rad
# endif
# ifdef UV_TIDES
                  UV_Tangle(i,j,itide)=UV_Tangle(i,j,itide)*deg2rad
                  UV_Tphase(i,j,itide)=UV_Tphase(i,j,itide)*deg2rad
# endif
                enddo
              enddo
              ierr=ierr_count
            else
              write(*,'(1x,4A)')  'WARNING: get_tides :: cannot find ',
     &       'variable ''tidal_period'' in ''', frcname(1:lfrc), '''.'
              if (ifrc < max_frc) then
                ifrc=ifrc+1
                goto 100
              endif
            endif
          enddo    !<-- itide
        endif
      endif  !<-- iic==ntstart
      return
      end
#else
      subroutine get_tides_empty
      end
#endif /* SSH_TIDES || UV_TIDES */
