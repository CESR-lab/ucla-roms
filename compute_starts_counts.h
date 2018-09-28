! Auxiliary module "compute_starts_counts.h":
!---------- ------ --------------------------
! Given two input variables, "type" and "record", which specify
! grid type of the field to be read/written from/to netCDF file, and
! record number (if there is an unlimited dimension for that field in
! the netCDF file), compute "start" and "count" arrays which contain
! starting indices and proper counts corresponding to the field data
! stored in the file.

! This module supports three different policies of reading/writing
! the data: (1) in the case of shared memory code the data is always
! written by single CPU into a single file, so that start/count are
! basically determined by the dimensions of arrays as it is written
! in the file;  (2) in the case of MPI code the data can still be
! written into a single file, or (3) multiple files containing data
! individually for each MPI-node (this option is activated by
! CPP-switch PARALLEL_FILES); In the case of MPI-single-file mode
! each MPI-node (one at-a-time) writes a rectangular block of array,
! so that "starts" depend on the position of the MPI-node on the
! processor grid, while "counts" corresponds to the size of MPI
! subdomains. In the case of PARALLEL_FILES "starts" corresponding
! to horizontal dimensions are both equal to 1, "counts" correspond
! to subdomains with physical boundaries kept, but MPI-halos striped
! out.

! Additionally, this module computes ranges "imin,imax,jmin,jmax"
! which define starting/ending indices of the portion of model array
! written into the file [hence in all cases count(1)=imax-imin+1, and
! count(2)=jmax-jmin+1].  This is necessary because model arrays
! have slightly different shapes than the corresponding arrays in the
! netCDF files.  The differences are due to index shifting (netCDF
! array index must always start from 1, while Fortran does not);
! Fortran array dimension padding; and stripping
! periodic/computational margins.

      integer imin,imax,jmin,jmax, start(4),count(4)

      jmin=horiz_type/2             ! calculate starting indices
      imin=horiz_type-2*jmin        ! in horizontal directions.

      ierr=0            ! These are default settings. In all cases
      do i=1,4          ! start,count(1:2) correspond to XI- and ETA-
        start(i)=1      ! dimensions, while 3 is either for vertical
        count(i)=1      ! dimension (if any) or for time record
      enddo             ! (2D-fields); 4 is for time record only

#ifdef MPI
      if (WESTERN_MPI_EDGE) then
        imin=imin + iwest-1       !<-- to account for grid type
      else
# ifndef PARALLEL_FILES
        start(1)=1+iSW_corn + 1-imin
# endif
        imin=iwest
      endif
      if (EASTERN_MPI_EDGE) then
        imax=ieast+1
      else
        imax=ieast
      endif
      if (SOUTHERN_MPI_EDGE) then
        jmin=jmin + jsouth-1
      else
# ifndef PARALLEL_FILES
        start(2)=1+jSW_corn + 1-jmin
# endif
        jmin=jsouth
      endif
      if (NORTHERN_MPI_EDGE) then
        jmax=jnorth+1
      else
        jmax=jnorth
      endif
#else                        /* non-MPI --> */
      imax=Lm+1 ; jmax=Mm+1
#endif
      count(1)=imax-imin+1
      count(2)=jmax-jmin+1

#ifndef NO_RECORD_CHECK
      if (nmax > 1) then
        count(3)=nmax
        start(4)=record
      else
        start(3)=record
      endif
#endif
