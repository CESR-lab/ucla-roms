! This is include file "tracers_defs.h"
!------ --- ---------------------------

! Passive tracers are defined here.

! iTandS represents the last index of the combination
! of 'temp' and 'salt'.
! Number passive tracers sequentially from iTandS+1

! Note: user must also put in index storing integer in
! the top of the tracers.h module for each passive tracer.
! E.g. if you add ptrace1, then at the top of tracers.h
! in the user input section, you need to add: integer iptrace1=0.

! (The reason for seperation is that the index variable needs to be
!  defined at compile time, and is not stored in an array like the
!  the rest of the variables that belong to each tracer.)

! itrace1 = integer to identify tracer in tracer array 't'
! wrt_t   = True/False whether to print tracer to output file
! t_vname = Tracer short name
! t_units = Tracer units (for outputing tracer)
! t_lname = Tracer long name (for outputing tracer)


      itrace1=1+iTandS;              wrt_t(itrace1) =.True.;
      t_vname(itrace1)='trace1';     t_units(itrace1)='%/%/%';
      t_tname(itrace1)='trace1_time'
      t_lname(itrace1)='long trace1'

      isalt2=2+iTandS;               wrt_t(isalt2) =.True.
      t_vname(isalt2)='salt2';       t_units(isalt2)='PSUuu'
      t_tname(isalt2)='salt2_time'
      t_lname(isalt2)='long salt2'

!      iptrace2=2+iTandS;            wrt_t(iptrace2) =.True.
!      t_vname(iptrace2)='ptrace2';  t_units(iptrace2)='uuunits'
!      t_tname(iptrace2)='ptrace2_time'
!      t_lname(iptrace2)='long ptrace2'
