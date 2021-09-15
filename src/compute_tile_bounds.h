! Auxiliary module "compute_tile_bounds.h":
!---------- ------ ------------------------
! Compute bounds designed to cover interior points of an array
! for shared memory subdomain partitioning (tiling.)
!
! input: tile -- usually from 0 to NSUB_X*NSUB_E-1 -- indicates
!                the specified subdomain;  tile=NSUB_X*NSUB_E
!                corresponds to the whole domain of RHO points
!                treated as a single block.
! outputs: istr,iend -- starting and ending indices of subdomain
!          jstr,jend    tile in XI- and ETA-directions.

      integer istr,iend, jstr,jend, i_X,j_E
#ifdef MPI
      integer size_X, margin_X, size_E, margin_E
#else
      integer, parameter :: size_X=(Lm+NSUB_X-1)/NSUB_X,
     &                    margin_X=(NSUB_X*size_X-Lm)/2,
     &                      size_E=(Mm+NSUB_E-1)/NSUB_E,
     &                    margin_E=(NSUB_E*size_E-Mm)/2
#endif

! 2021/03: currently ALLOW_SINGLE_BLOCK_MODE is not defined.
#ifdef ALLOW_SINGLE_BLOCK_MODE
C$    integer trd, omp_get_thread_num
      if (tile==NSUB_X*NSUB_E) then
C$      trd=omp_get_thread_num()
C$      if (trd>0) return !--> just return, if not master thread
# ifdef MPI
        istr=iwest      ! MONOBLOCK VERSION:
        iend=ieast      ! Do not divide grid
        jstr=jsouth     ! into tiles.
        jend=jnorth
# else
        istr=1
        iend=Lm
        jstr=1
        jend=Mm
      else
# endif
#endif

        j_E=tile/NSUB_X
        i_X=tile-j_E*NSUB_X
        if (mod(j_E,2)==1) i_X=NSUB_X-1 -i_X   !<-- sweep reversal

#ifdef MPI
        if (mod(inode,2)>0) then               ! make mirror-image
          i_X=NSUB_X-1 -i_X                    ! symmetry for sweep
        endif                                  ! trajectories of MPI
        if (mod(jnode,2)>0) then               ! subdomains adjacent
          j_E=NSUB_E-1 -j_E                    ! in both directions
        endif

        size_X=(ieast-iwest+NSUB_X)/NSUB_X
        margin_X=(NSUB_X*size_X - ieast+iwest-1)/2
        istr=iwest-margin_X + i_X*size_X
        iend=min( istr + size_X-1 ,ieast)
        istr=max(istr,iwest)
#else
        istr=1-margin_X + i_X*size_X
        iend=min( istr + size_X-1, Lm)
        istr=max(istr,1)
#endif


#ifdef MPI
        size_E=(jnorth-jsouth +NSUB_E)/NSUB_E
        margin_E=(NSUB_E*size_E -jnorth+jsouth-1)/2
        jstr=jsouth-margin_E + j_E*size_E
        jend=min( jstr + size_E-1 ,jnorth)
        jstr=max(jstr,jsouth)
#else
        jstr=1-margin_E + j_E*size_E
        jend=min( jstr + size_E-1, Mm)
        jstr=max(jstr,1)
#endif


#ifdef ALLOW_SINGLE_BLOCK_MODE
      endif
#endif

