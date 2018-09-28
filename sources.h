#ifdef PSOURCE

! Nsrc       Number of point Sources/Sinks.
! Dsrc       Direction of point Sources/Sinks:  0 = along XI-;
!                                          1= along ETA-direction.
! Isrc,Jsrc  i,j-grid location of point Sources/Sinks,
! Qsrc       Mass transport profile (m3/s) of point Sources/Sinks.
! Qbar       Vertically integrated Qsrc (m3/s) of point
! QbarG      Latest two-time snapshots of vertically integrated
!              mass transport (m3/s) of point Sources/Sinks.
! Tsrc       Tracer (tracer units) point Sources/Sinks.
! TsrcG      Latest two-time snapshots of tracer (tracer units)
!              point Sources/Sinks.
! Qshape     Nondimensional shape function to distribute mass point
!            Sources/Sinks vertically.

      real Qbar(Msrc), Qsrc(Msrc,N), Qshape(Msrc,N), Tsrc(Msrc,N,NT)
      integer Nsrc, Dsrc(Msrc), Isrc(Msrc), Jsrc(Msrc)
      common /psources/Qbar,Qsrc,Qshape, Tsrc, Nsrc, Dsrc,Isrc,Jsrc


# ifndef ANA_PSOURCE
#  if defined PSOURCE_DATA || defined ALL_DATA
#   undef PSOURCE_DATA
      real QbarG(Msrc,2), TsrcG(Msrc,2,NT)
      real(kind=8) psrc_cycle, psrc_time(2)
      integer psrc_ncycle, psrc_rec, itpsrc, ntpsrc,
     &        ncidpsrcs, psrc_tid,Qbar_id,Tsrc_id(NT)
      common /psrcs_data/
     &        QbarG, TsrcG, psrc_cycle, psrc_time,
     &        psrc_ncycle, psrc_rec, itpsrc, ntpsrc,
     &        ncidpsrcs, psrc_tid,Qbar_id,Tsrc_id
      character(len=64) psrc_file
      common /psrcs_data_file/ psrc_file
#  endif
# endif
#endif
