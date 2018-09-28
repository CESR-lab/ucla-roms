! Variables for internal use by message exhanging routines:
!---------- --- -------- --- -- ------- --------- ---------
! inode,jnode -- indices of current MPI node on processor grid;
! p_W,..,p_NW -- MPI ranks (within ocean_grid_comm) of the eight
!                neighboring subdomains;
! exc_call_count -- counter for mpi_exchange calls;
! west_msg_exch,...,north_msg_exch -- logical switch to send/receive
!                MPI messages to/from the corresponding directions
 
#ifdef MPI
      integer inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
#endif
 
