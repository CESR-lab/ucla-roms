#include "cppdefs.opt"
                                           ! Checking compatibility of
      subroutine check_tiling (ierr)       ! dimensions of the grid and
      use param
      implicit none                        ! partitioning parameters:
      integer ierr, tile                   ! logic of the code requires

      do tile=0,NSUB_X*NSUB_E-1            ! all the tiles to be at
        call check_tile_bounds (tile,ierr) ! least two point wide.
      enddo                                ! Setting NSUB_X and NSUB_E
      end                                  ! incorrectly may result in
                                           ! smaller size.

      subroutine check_tile_bounds (tile, ierr)
      use param
      use hidden_mpi_vars
      implicit none
      integer tile, ierr, isize,jsize
#include "compute_tile_bounds.h"

      isize=iend-istr+1                    ! NOTE: run terminates, if
      jsize=jend-jstr+1                    ! violation occurs here.

      if (isize<2) then
        ierr=ierr+1
        write(*,'(1x,2A,I3,1x,A,I4,3x,A,I3,1x,A,2I4)') '### ERROR: ',
     &          'NSUB_X =', NSUB_X,  'is incompatible with Lm =', Lm,
     &          'i_X =', i_X,  '==> istr,iend =', istr,iend
      endif
      if (jsize<2) then
        ierr=ierr+1
        write(*,'(1x,2A,I3,1x,A,I4,3x,A,I3,1x,A,2I4)') '### ERROR: ',
     &          'NSUB_E =', NSUB_E,  'is incompatible with Mm =', Mm,
     &          'j_E =', j_E,  '==> jstr,jend =', jstr,jend
      endif
      end

      subroutine print_switches

! Analyze content of string "cpps" and print out all filenames
! (substrings enclosed in <..>) and CPP switches activated in these
! files as they were recorded in "cpps" by several other routines,
! starting with "check_switches1".    This routine is merely for
! reporting, and it does not perform any checks, nor does it raise
! the error flag.
      use strings
      implicit none
      character(len=32) bffr(3)
      integer is,ie, i,nwd, lbffr(3)

      write(*,'(/1x,A)') 'Activated CPP configuration switches:'
      nwd=0
      is=1
      do while (is<max_opt_size)
        do while (is<max_opt_size .and. cpps(is:is)==' ')
          is=is+1
        enddo
        ie=is
        do while (ie<max_opt_size .and. cpps(ie:ie)/=' ')
          ie=ie+1
        enddo
        if (cpps(ie:ie) == ' ') ie=ie-1
        if (cpps(is:is) == '<' .and. cpps(ie:ie) == '>') then
          if (nwd>0) then
            write(*,1) (bffr(i)(1:lbffr(i)), i=1,nwd)
            nwd=0
          endif
          write(*,'(1x,A)') cpps(is:ie)
        else
          nwd=nwd+1
          bffr(nwd)=cpps(is:ie)
          lbffr(nwd)=ie-is+1
          if (nwd==3) then
            write(*,1) (bffr(i)(1:lbffr(i)), i=1,nwd)
            nwd=0
          endif
        endif
        is=ie+1                      ! This format is designed to
      enddo                          ! accommodate three columns at
  1   format(T8,A,T32,A,T56,A)       ! most 24-character-wide and fit
      end                            ! it into 80-character width.


      subroutine check_switches2(ierr)

! Manually check mutual consistency of settings of CPP switches, report
! contradictions (if any), and raise errors flag to terminate the run.

      use param
      use scalars
      use strings
      implicit none
      integer ierr, nobc_west, nobc_east, nobc_south, nobc_north,
     &                                ncoupl, nAk, is,ie, lenstr

      ie=lenstr(cpps)

! Track settings in "set_global_definitions.h". The following
! section is generated automatically by executing
!
!            cppcheck set_global_definitions.h
!
! and copy-pasting a manually selected portion of "check_switches1.F"
! generated this way. This is preferred over just letting cppcheck to
! automatically generate the code for tracking switches defined in
! "set_global_definitions.h" because most of what the latter contains
! is a set of permanently defined macros tracking which which would
! only clutter the signature stored in "cpps" with adding little
! useful information.

!----->>>>>>>>> begin copy-pasted section.

      is=ie+2 ; ie=is+25
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='<set_global_definitions.h>'
      cpps(ie+1:ie+1)=' '
#ifdef PRED_COUPLED_MODE
      is=ie+2 ; ie=is+16
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='PRED_COUPLED_MODE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CORR_COUPLED_MODE
      is=ie+2 ; ie=is+16
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CORR_COUPLED_MODE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EXTRAP_BAR_FLUXES
      is=ie+2 ; ie=is+16
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EXTRAP_BAR_FLUXES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef KEEP_CORIOLIS
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='KEEP_CORIOLIS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef IMPLCT_NO_SLIP_BTTM_BC
      is=ie+2 ; ie=is+21
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='IMPLCT_NO_SLIP_BTTM_BC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VAR_RHO_2D
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='VAR_RHO_2D'
      cpps(ie+1:ie+1)=' '
#endif

! <<<<<<<<-- end of copy-pasted section.

! Checking correctness of coupling flags:

      ncoupl=0
#ifdef PRED_COUPLED_MODE
      ncoupl=ncoupl+1
#endif
#ifdef CORR_COUPLED_MODE
      ncoupl=ncoupl+1
#endif
      if (ncoupl==0) then
        write(*,'(/1x,3A/)') '### ERROR: neither PRED_COUPLED_MODE ',
     &                       'nor CORR_COUPLED_MODE are defined in ',
     &                       '"set_global_definitions.h".'
        ierr=ierr+1
      elseif (ncoupl>1) then
        write(*,'(/1x,2A/12x,2A/)')     '### ERROR: miss-configured ',
     &          'switches PRED_COUPLED_MODE and CORR_COUPLED_MODE ',
     &          'in "set_global_definitions.h": only one of them ',
     &          'should be defined.'
        ierr=ierr+1
      endif

! Checking the number of selections of boundary conditions.

#ifdef EW_PERIODIC
      nobc_west=1     ! here "nobc_west=1" means that it is already
      nobc_east=1     ! chosen and no further selection is possible
#else
      nobc_west=0
      nobc_east=0
#endif
#ifdef OBC_WEST
      nobc_west=nobc_west+1
#endif
#ifdef OBC_EAST
      nobc_east=nobc_east+1
#endif
      if (nobc_west>1) then
        write(*,'(1x,2A)') '### ERROR in "cppdefs.opt": more than one ',
     &          'boundary condition is chosen on the WESTERN EDGE.'
        ierr=ierr+1
      endif
      if (nobc_east>1) then
        write(*,'(1x,2A)') '### ERROR in "cppdefs.opt": more than one ',
     &          'boundary condition is chosen on the EASTERN EDGE.'
        ierr=ierr+1
      endif

#ifdef NS_PERIODIC
      nobc_south=1
      nobc_north=1
#else
      nobc_south=0
      nobc_north=0
#endif
#ifdef OBC_SOUTH
      nobc_south=nobc_south+1
#endif
#ifdef OBC_NORTH
      nobc_north=nobc_north+1
#endif
      if (nobc_south>1) then
        write(*,'(1x,2A)') '### ERROR in "cppdefs.opt": more than one ',
     &          'boundary condition is chosen on the SOUTHERN EDGE.'
        ierr=ierr+1
      endif
      if (nobc_north>1) then
        write(*,'(1x,2A)') '### ERROR in "cppdefs.opt": more than one ',
     &          'boundary condition is chosen on the NORTHERN EDGE.'
        ierr=ierr+1
      endif

! Checking the number of selections vertical mixing schemes.

      nAk=0
#ifdef ANA_VMIX
      nAk=nAk+1
#endif
#ifdef BVF_MIXING
      nAk=nAk+1
#endif
#ifdef LMD_MIXING
      nAk=nAk+1
#endif
#ifdef MY2_MIXING
      nAk=nAk+1
#endif
#ifdef MY25_MIXING
      nAk=nAk+1
#endif
#ifdef PP_MIXING
      nAk=nAk+1
#endif
      if (nAk>1) then
        write(*,'(1x,2A)') '### ERROR in "cppdefs.opt": more than ',
     &                     'one vertical mixing scheme is chosen.'
        ierr=ierr+1
      endif

/*
#ifdef SPLIT_EOS
# ifndef NONLIN_EOS
      write(*,'(1x,2A/12x,A)') '### ERROR in "cppdefs.opt": ',
     &        'illegal combination of CPP-switches: ',
     &        'SPLIT_EOS is defined, but NONLIN_EOS is not.'
      ierr=ierr+1
# endif
#endif
*/

! Signal to quit if explicit time-step splitting on shallow water
! set-up.

#ifndef SOLVE3D
      if (ndtfast>1) then
        write(*,'(1x,A,I3,1x,A/12x,A,2x,A)')  '### ERROR: NDTFAST =',
     &         ndtfast, 'is greater than unity for a shallow water',
     &        'configuration.','Change it to unity in startup file.'
     &
        ierr=ierr+1
      endif
#endif
      return
  99  mpi_master_only write(*,'(/1x,2A/14x,A)')
     &  '### ERROR: checkdefs :: Insufficient size of string "cpps" ',
     &  'in file "strings".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
