mg_autotune.o : mg_autotune.f90 mg_solvers.o mg_projection.o mg_tictoc.o mg_namelist.o mg_mpi.o mg_cst.o
mg_cst.o : mg_cst.f90
mg_gather.o : mg_gather.f90 mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_grids.o : mg_grids.f90 mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_horiz_grids.o : mg_horiz_grids.f90 mg_namelist.o mg_gather.o mg_mpi_exchange.o mg_grids.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_intergrids.o : mg_intergrids.f90 mg_gather.o mg_mpi_exchange.o mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_mpi_exchange.o : mg_mpi_exchange.f90 mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_mpi.o : mg_mpi.f90 mg_cst.o
mg_namelist.o : mg_namelist.f90 mg_tictoc.o mg_cst.o
mg_netcdf_out.o : mg_netcdf_out.f90 mg_cst.o
mg_projection.o : mg_projection.f90 mg_netcdf_out.o mg_gather.o mg_mpi_exchange.o mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_relax.o : mg_relax.f90 mg_mpi_exchange.o mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_solvers.o : mg_solvers.f90 mg_netcdf_out.o mg_relax.o mg_intergrids.o mg_grids.o mg_namelist.o mg_tictoc.o mg_mpi.o mg_cst.o
mg_tictoc.o : mg_tictoc.f90
mg_vert_grids.o : mg_vert_grids.f90 mg_netcdf_out.o mg_namelist.o mg_gather.o mg_mpi_exchange.o mg_grids.o mg_tictoc.o mg_mpi.o mg_cst.o
nhmg_debug.o : nhmg_debug.f90 mg_netcdf_out.o mg_mpi.o
nhmg.o : nhmg.f90 mg_netcdf_out.o mg_solvers.o mg_projection.o mg_vert_grids.o mg_horiz_grids.o mg_autotune.o mg_mpi_exchange.o mg_tictoc.o mg_namelist.o mg_grids.o mg_mpi.o mg_cst.o
