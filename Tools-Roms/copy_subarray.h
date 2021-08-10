      integer vdims, start(vdims), count(vdims), count1(vdims)

! This file is needed by "partit.F" and "ncjoin.F" and it contains
! declaration of dimensional variables packed together with a segment
! of executable code to copy rectangular block of data "buff" (which is
! just an MPI subdomain) into array "bfr_out" which covers the entire
! physical grid, or vice versa, extracts "buff" from the whole-grid
! array "bfr_in" into "buff".  Arrays "buff","bfr_in", "bfr_out"
! themselves are not declared here because they may be of different
! type: integer of kind=1(byte), 2, or 4, real(kind=4 or 8), or
! character(len=1), so they are rather declared directly in "partit.F"
! or "ncjoin.F" where this file is used to create multiple versions of
! alike subroutines to process different kinds of input data.
! Essentially extraction of this code segment into a separate file
! is poor-man polymorphism. 

! The input arguments (declared above in this file) are:

! vdims -- the actual number of logical dimensions of arrays "buff",
!                                            "bfr_in", and "bfr_out';

! start -- an array containing starting indices of the block of data
!                          in index-space of the whole physical grid;

! count -- an array containing sizes of the block, essentially as it
!          would be dimensioned as

!               buff(count(1),count(2),count(3), etc...)

!          except that formally speaking "buff" is declared as
!          1D-array because the actual number of its dimensions
!          "vdims" may be anything during run time; 

! count1 -- an array containing dimensions of the entire grid, as it
!           would be

!              bfr_in,bfr_out(count1(1),count1(2),count1(3), etc...)


      integer ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1

! Check whether the last dimension is unlimited dimension and if so,
! ignore it, since records corresponding to the different indices
! along the unlimited dimension are always written one-by-one.

      if (count1(vdims)==1) then
        ndims=vdims-1             ! WARNING: This code is restricted
      else                        ! to have no more than 4 dimensions
        ndims=vdims               ! for partitioned arrays not counting
      endif                       ! unlimited dimension.
      istr=start(1)
      imax=count(1)               ! Furthermore, it is assumed that
      imax1=count1(1)             ! unlimited dimension is always the
      if (ndims > 1) then         ! last one [this is the standard
        jstr=start(2)             ! practice in ROMS, however netCDF
        jmax=count(2)             ! is not restricted to do so].
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif

c*    write(*,'(1x,A,2I3,3(3x,A,4I4))')  'ndims,vdims =', ndims,
c*   &   vdims, 'imax1,jmax1,kmax1,lmax1 =', imax1,jmax1,kmax1,lmax1,
c*   &                  'imax,jmax,kmax,lmax =', imax,jmax,kmax,lmax,
c*   &                  'istr,jstr,kstr,lstr =', istr,jstr,kstr,lstr

      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
#if defined PARTIT
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
#elif defined NCJOIN
              bfr_out(i+istr-1 +js1*imax1) = buff(i +js*imax)
#else
              ??????
#endif
            enddo
          enddo
        enddo
      enddo
