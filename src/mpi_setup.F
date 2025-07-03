#include "cppdefs.opt"
#ifdef MPI
c--# define TRANSPOSED_MPI_NODE_ORDER
c--# define PRINT_HOSTNAME

      subroutine mpi_setup(ierr)
      use param
      use hidden_mpi_vars
      use mpi
      use roms_read_write
      use netcdf 
      use nc_read_write
      implicit none
      integer :: ierr,newpart,ncid
      integer,dimension(8) :: neighbors
      integer,dimension(2) :: subdompos
      integer,dimension(4) :: partition 
      character(len=128) :: fname

      integer off_xi,off_eta, conv_ijnode_to_rank,check_neighbor
# ifdef PRINT_HOSTNAME
      integer status(MPI_STATUS_SIZE), blank
      character(len=38) msg1 ! <-- greeting message to print
      character(len=23) msg2 
# endif

      exc_call_count=0

      ! read in gridfile to check if mpi masking is present
      ierr=nf90_open(grdname,nf90_nowrite,ncid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr)
      newpart=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
      ! determine nodes based on version
      if (newpart.eq.nf90_noerr) then
        ierr=nf90_get_att(ncid,nf90_global,"partition",partition)
        NNODES=partition(2)
        ! grab neighbors & i/j nodes
        ierr=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
        ierr=nf90_get_att(ncid,nf90_global,"subdompos",subdompos)
        print*,'MPI_SETUP: MASKING ON'
      else
        NNODES=NP_XI*NP_ETA
      endif
      ierr=nf90_close(ncid)


      if (nsize==NNODES) then       ! Check whether the number of
                                    ! nodes specified by -np argument
# ifdef TRANSPOSED_MPI_NODE_ORDER
        inode=mynode/NP_ETA         ! of mpiexec matches the parameter
        jnode=mynode-inode*NP_ETA   ! settings in the code, and find
# else
        if (newpart.eq.nf90_noerr) then
        inode=subdompos(1)
        jnode=subdompos(2)
        else
        jnode=mynode/NP_XI          ! indices inode,jnode identifying
        inode=mynode-jnode*NP_XI    ! the location of current subdomain
        endif
# endif
                                    ! on the "processor grid".
# ifdef EW_PERIODIC
        west_msg_exch=.true.        ! Determine whether the subdomain
        east_msg_exch=.true.        ! has neighbors on the four sides
# else
        if (newpart.eq.nf90_noerr) then
        west_msg_exch=(neighbors(7)/=-99)
        east_msg_exch=(neighbors(3)/=-99)
        else
        west_msg_exch=inode>0            ! around it and set the
        east_msg_exch=inode<NP_XI-1      ! corresponding logical
        endif
# endif

# ifndef EW_PERIODIC
        west_exchng=west_msg_exch   ! receive incoming  messages
        east_exchng=east_msg_exch   ! to/from the neighbor.
# endif

# ifdef NS_PERIODIC
        south_msg_exch=.true.       ! via exchange of computational
        north_msg_exch=.true.       ! margins, so that in this case
# else
        if (newpart.eq.nf90_noerr) then
        south_msg_exch=(neighbors(5)/=-99)
        north_msg_exch=(neighbors(1)/=-99)
        else
        south_msg_exch=jnode>0         ! communications take place
        north_msg_exch=jnode<NP_ETA-1  ! even if the subdomain is
        endif
# endif

# ifndef NS_PERIODIC
        south_exchng=south_msg_exch
        north_exchng=north_msg_exch
# endif

! Determine MPI-ranks of the MPI subdomains adjacent from the sides
! and corners. These are needed to designate MPI ranks for sources of
! incoming and targets for outgoing messages.

        if (newpart.eq.nf90_noerr) then
         p_E=check_neighbor(neighbors(3))
        p_NE=check_neighbor(neighbors(2))
        p_N =check_neighbor(neighbors(1))
        p_NW=check_neighbor(neighbors(8))
         p_W=check_neighbor(neighbors(7))
        p_SW=check_neighbor(neighbors(6))
        p_S =check_neighbor(neighbors(5))
        p_SE=check_neighbor(neighbors(4))
        else
         p_E=conv_ijnode_to_rank(inode+1,jnode  , NP_XI,NP_ETA)
        p_NE=conv_ijnode_to_rank(inode+1,jnode+1, NP_XI,NP_ETA)
        p_N =conv_ijnode_to_rank(inode  ,jnode+1, NP_XI,NP_ETA)
        p_NW=conv_ijnode_to_rank(inode-1,jnode+1, NP_XI,NP_ETA)
         p_W=conv_ijnode_to_rank(inode-1,jnode  , NP_XI,NP_ETA)
        p_SW=conv_ijnode_to_rank(inode-1,jnode-1, NP_XI,NP_ETA)
        p_S =conv_ijnode_to_rank(inode  ,jnode-1, NP_XI,NP_ETA)
        p_SE=conv_ijnode_to_rank(inode+1,jnode-1, NP_XI,NP_ETA)
        endif

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

! Convert masked value, if any, to -1.

      function check_neighbor(n)
      implicit none
      integer check_neighbor,n
      if (n==-99) then
        check_neighbor=-1
      else
        check_neighbor=n
      endif
      end 

#else
      subroutine mpi_setup_empty
      end
#endif    /* MPI */
