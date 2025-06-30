#include "cppdefs.opt"
#ifdef MPI
c--# define TRANSPOSED_MPI_NODE_ORDER
c--# define PRINT_HOSTNAME

C$    subroutine master_num_threads (nthrds)
C$    implicit none
C$    integer nthrds, trd, omp_get_num_threads, omp_get_thread_num
C$    trd=omp_get_thread_num()
C$    if (trd == 0) nthrds=omp_get_num_threads()
C$    end


      subroutine mpi_setup(ierr)
      use param
      use hidden_mpi_vars
      use mpi
      use roms_read_write
      implicit none
      integer ierr

C$    integer nthrds
      integer nsize, off_xi,off_eta, conv_ijnode_to_rank
# ifdef PRINT_HOSTNAME
      integer status(MPI_STATUS_SIZE), blank
      character(len=38) msg1 ! <-- greeting message to print
      character(len=23) msg2
# endif

      ocean_grid_comm=MPI_COMM_WORLD
      call MPI_Comm_size(ocean_grid_comm, nsize,  ierr)
      call MPI_Comm_rank(ocean_grid_comm, mynode, ierr)

      exc_call_count=0
                                         ! Determine number of threads
C$    if (mynode == 0) then              ! on MPI node rank=0, then
C$OMP PARALLEL SHARED(nthrds)            ! broadcast it, so all others
C$      call master_num_threads (nthrds) ! can set the same number of
C$OMP END PARALLEL                       ! threads as the rank=0 node.
C$    endif
C$    call MPI_Bcast(nthrds, 1, MPI_INTEGER, 0, ocean_grid_comm, ierr)
C$    if (mynode > 0) then
C$      call omp_set_num_threads (nthrds)
C$    endif
                                    ! Check whether the number of
      if (nsize == NNODES) then     ! nodes specified by -np argument
# ifdef TRANSPOSED_MPI_NODE_ORDER
        inode=mynode/NP_ETA         ! of mpiexec matches the parameter
        jnode=mynode-inode*NP_ETA   ! settings in the code, and find
# else
        jnode=mynode/NP_XI          ! indices inode,jnode identifying
        inode=mynode-jnode*NP_XI    ! the location of current subdomain
# endif
                                    ! on the "processor grid".
!       if (NP_XI > 1) then
# ifdef EW_PERIODIC
          west_msg_exch=.true.      ! Determine whether the subdomain
          east_msg_exch=.true.      ! has neighbors on the four sides
# else
          west_msg_exch=inode>0          ! around it and set the
          east_msg_exch=inode<NP_XI-1    ! corresponding logical
# endif
!       else                        ! flags: e.g. est_msg_exch==.true.
!         west_msg_exch=.false.     ! means that there is a neighbor on
!         east_msg_exch=.false.     ! the west side, so this subdomain
!       endif                       ! must send and and expects to
# ifndef EW_PERIODIC
        west_exchng=west_msg_exch   ! receive incoming  messages
        east_exchng=east_msg_exch   ! to/from the neighbor.
# endif
                                    ! Note that periodic boundaries
!       if (NP_ETA > 1) then        ! (if any) are treated exclusively
# ifdef NS_PERIODIC
          south_msg_exch=.true.     ! via exchange of computational
          north_msg_exch=.true.     ! margins, so that in this case
# else
          south_msg_exch=jnode>0         ! communications take place
          north_msg_exch=jnode<NP_ETA-1  ! even if the subdomain is
# endif
!       else                        ! located near the edge of
!         south_msg_exch=.false.    ! the grid.
!         north_msg_exch=.false.
!       endif
# ifndef NS_PERIODIC
        south_exchng=south_msg_exch
        north_exchng=north_msg_exch
# endif

! Determine MPI-ranks of the MPI subdomains adjacent from the sides
! and corners. These are needed to designate MPI ranks for sources of
! incoming and targets for outgoing messages.

         p_E=conv_ijnode_to_rank(inode+1,jnode  , NP_XI,NP_ETA)
        p_NE=conv_ijnode_to_rank(inode+1,jnode+1, NP_XI,NP_ETA)
        p_N =conv_ijnode_to_rank(inode  ,jnode+1, NP_XI,NP_ETA)
        p_NW=conv_ijnode_to_rank(inode-1,jnode+1, NP_XI,NP_ETA)
         p_W=conv_ijnode_to_rank(inode-1,jnode  , NP_XI,NP_ETA)
        p_SW=conv_ijnode_to_rank(inode-1,jnode-1, NP_XI,NP_ETA)
        p_S =conv_ijnode_to_rank(inode  ,jnode-1, NP_XI,NP_ETA)
        p_SE=conv_ijnode_to_rank(inode+1,jnode-1, NP_XI,NP_ETA)

!       print *, 'neighbors: ', mynode, p_W,p_E

! Determine bounds of the usable portion of model arrays:

        off_xi=NP_XI*Lm-LLm
        if (inode == 0) then               ! DevinD added
          iSW_corn = 0                     ! need to fix iSW_corn as used to determine WESTERN_EDGE
        else
          iSW_corn = inode*Lm - off_xi/2   ! add length of western domain.
        endif
        if (inode == 0) then
          Lm = Lm - off_xi/2
        endif
        if (inode == NP_XI-1) then
          Lm = Lm -(off_xi+1)/2
        endif
        iwest = 1
        ieast = Lm                               ! all sub-domains (Lm same for all interior,
                                                 ! possibly different E & W edges

        off_eta=NP_ETA*Mm-MMm
        jSW_corn=jnode*Mm-off_eta/2
!        if (jnode == 0) then
!          jsouth=1+off_eta/2
!        else
!          jsouth=1
!        endif
!        if (jnode < NP_ETA-1) then
!          jnorth=Mm
!        else
!          jnorth=Mm -(off_eta+1)/2
!        endif
        if (jnode == 0) then                ! DevinD added
          jSW_corn = 0                      ! need to fix iSW_corn as used to determine WESTERN_EDGE
        else
          jSW_corn = jnode*Mm - off_eta/2   ! add length of western domain.
        endif
        if (jnode == 0) then
          Mm = Mm - off_eta/2
        endif
        if (jnode == NP_ETA-1) then
          Mm = Mm -(off_eta+1)/2
        endif
        jsouth = 1
        jnorth = Mm

#ifdef VERBOSE
      write(*,'(A,7I5,1x,A,I4)') 'XI:', LLm, off_xi, iSW_corn, Lm,
     & ieast-iwest+1, iwest+iSW_corn,ieast+iSW_corn, 'node=', mynode
      write(*,'(A,7I5,1x,A,I4)') 'ETA:',MMm, off_eta, jSW_corn, Mm,
     & jnorth-jsouth+1,jsouth+jSW_corn,jnorth+jSW_corn,'node=',mynode
#endif

# ifdef PARALLEL_FILES
         xi_rho=ieast-iwest+1
         if (EASTERN_MPI_EDGE) then
           xi_rho=xi_rho+1
         endif
         if (WESTERN_MPI_EDGE) then
           xi_rho=xi_rho+1
           xi_u=xi_rho-1
         else
           xi_u=xi_rho
         endif

         eta_rho=jnorth-jsouth+1
         if (NORTHERN_MPI_EDGE) then
           eta_rho=eta_rho+1
         endif
         if (SOUTHERN_MPI_EDGE) then
           eta_rho=eta_rho+1
           eta_v=eta_rho-1
         else
           eta_v=eta_rho
         endif
# endif
! Print hostnames for each rank = 0 to NNODES-1 in the sequence.
# ifdef PRINT_HOSTNAME
        call MPI_Barrier(ocean_grid_comm, ierr)
        if (mynode>0) call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,1,
     &                                   ocean_grid_comm, status, ierr)

        write(msg1,'(A,I5,1x,A,I5,1x,A,2x)') 'MPI rank', mynode,
     &                              'out of', nsize, 'running on'
        write(msg2,'(2x,A,2I4)') 'inode,jnode =', inode,jnode

        call system('echo "'/ /msg1/ /'`hostname -s`'/ /msg2/ /'"')
        if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                          mynode+1, 1, ocean_grid_comm, ierr)
        call MPI_Barrier(ocean_grid_comm, ierr)
# endif

        ierr=0
      else
        mpi_master_only write(*,'(/1x,A,I4,1x,A,I3,A/)')
     &   '### ERROR: mpi_setup: number of MPI-nodes should be',
     &                         NNODES, 'instead of ', nsize, '.'
        ierr=99
      endif
      end

! Convert processor-grid indices i,jnode into MPI rank. Check whether
! the indices are within the bounds of processor grid. Return them back
! into periodicity domain, if there is periodicity in either direction,
! but only if the number of subdomains is greater than one (hence
! periodicity must go through message passing).

      function conv_ijnode_to_rank(ii,jj, np_xi,np_eta)
      implicit none
      integer conv_ijnode_to_rank, ii,jj, np_xi,np_eta, i,j
      i=ii ; j=jj
# ifdef EW_PERIODIC
      if (i < 0) then
        i=i+np_xi
      elseif (i > np_xi-1) then
        i=i-np_xi
      endif
# endif
# ifdef NS_PERIODIC
      if (j < 0) then
        j=j+np_eta
      elseif (j > np_eta-1) then
        j=j-np_eta
      endif
# endif
      if (0 <= i .and. i < np_xi .and. 0 <= j .and. j < np_eta) then
# ifdef TRANSPOSED_MPI_NODE_ORDER
        conv_ijnode_to_rank=j + i*np_eta
# else
        conv_ijnode_to_rank=i + j*np_xi
# endif
      else
        conv_ijnode_to_rank=-1 !<- meaning that the rank does not exist
      endif
      end

#else
      subroutine mpi_setup_empty
      end
#endif    /* MPI */
