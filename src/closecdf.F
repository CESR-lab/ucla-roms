#include "cppdefs.opt"

      subroutine closecdf       ! Close output netCDF files.
      use param
      use scalars
      implicit none

# ifdef MPI_SILENT_MODE
      if (mynode == 0) then
# endif
        write(*,'(/1x,A/)') 'MAIN: DONE - Records written:'

        if (may_day_flag==1) then
          write(*,'(/1x,A/)') 'MAIN: Abnormal termination: BLOWUP'
        elseif (may_day_flag==2) then
          write(*,'(/1x,2A/)')  '### ERROR: Cause of termination: ',
     &                                            'netCDF input.'
        elseif (may_day_flag==3) then
          write(*,'(/1x,2A/)')  '### ERROR: Cause of termination: ',
     &                                           'netCDF output.'
        elseif (may_day_flag==4) then
          write(*,'(/1x,A/)') '### ERROR: Cannot open netCDF file.'
        elseif (may_day_flag==5) then
          write(*,'(/1x,A/)') '### ERROR: Error in input script.'
        elseif (may_day_flag==6) then
          write(*,'(/1x,2A/)') '### ERROR: read_inp :: An input ',
     &                                   'file is not available.'
        elseif (may_day_flag==7) then
          write(*,'(/1x,A/)') '### ERROR: KPP algoritm failure.'
        elseif (may_day_flag==8) then
          write(*,'(/1x,A/)') '### ERROR: Wrong choice of tiling.'
        elseif (may_day_flag==9) then
          write(*,'(/1x,2A/)') '### ERROR: Wrong choice of ',
     &                               'S-coordinate parameters.'
        endif
# ifdef MPI_SILENT_MODE
      endif
# endif
      end
