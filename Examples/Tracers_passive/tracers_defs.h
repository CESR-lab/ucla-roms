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

! itrace1   = integer to identify tracer in tracer array 't'
! wrt_t     = True/False whether to print tracer to output file
! t_vname   = Tracer short name
! t_units   = Tracer units (for outputing tracer)
! t_lname   = Tracer long name (for outputing tracer)
! t_ana_frc = Whether surf flux is read in (0), or analytical (1)
!             Could be extended 2,3,... depending on different types of
!             analytical forcing, e.g. time variant/invariant, and further
!             specific for each variable using its 'itrace' index...
! itot      = Counter to increment tracer indices without hardcoding indices,
!             thus each new tracer index (e.g. itrcX) simply needs:
!             itrcX=itot;    and   itot=itot+1 afterwards for the next tracer.

! wrt_t_diag= Only with DIAGNOSTICS flag. Don't need otherwise.
!             True/False to output tracer diagnostics.



! - TEMP & SALT:

!     Defined in tracers.F, but outputting controller here:
      wrt_t(itmp) =.True.;           wrt_t_avg(itmp) =.True.
#ifdef SALINITY
      wrt_t(islt) =.True.;           wrt_t_avg(islt) =.True.
#endif



! - PASSIVE TRACERS:

      itrace1=itot;                   itot=itot+1
      wrt_t(itrace1)  =.True.;        wrt_t_avg(itrace1)=.True.
      t_vname(itrace1)='trace1';      t_units(itrace1)  ='%/%/%'
      t_tname(itrace1)='trace1_time'; t_ana_frc(itrace1)=0
      t_lname(itrace1)='long trace1'

      isalt2=itot;                    itot=itot+1               
      wrt_t(isalt2) =.True.;          wrt_t_avg(isalt2)=.True.
      t_vname(isalt2)='salt2';        t_units(isalt2)  ='PSUuu'
      t_tname(isalt2)='salt2_time';   t_ana_frc(isalt2)=0
      t_lname(isalt2)='long salt2'



! - BGC TRACERS:



