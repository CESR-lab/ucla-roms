# Installation

## Installing dependencies
ROMS runs on Unix-like systems. It additionally requires:

- A Fortran compiler (`gfortran` and `ifx` are actively supported)
- The `netcdf-fortran` library for netCDF i/o
- An MPI library (e.g. MPICH, OpenMPI, Intel MPI) for parallelization

... Using ROMS with MARBL biogeochemistry additionally requires MARBL v0.45 (TODO)

### Obtaining dependencies via conda or containers:
⚠️ **NOTE** : this option is best-suited to small-scale scientific computing systems, cloud instances, and personal computers/laptops. Large-scale supercomputers typically have internal package management systems that vary between machines and over time. ⚠️

The necessary (GNU) dependencies for ROMS can be installed via conda in an environment created from an included `yaml` file:

```
conda env create -f ci/environment.yml --name roms_env
conda activate roms_env
```

TODO: set environment variables

Alternatively, ROMS' dependencies come pre-built in containers:

- [Image 1](https://ghcr.io/dafyddstephenson/roms_gfortran_build_env:1.0) (uses gfortran, MPICH)
- [Image 2](https://ghcr.io/cworthy-ocean/marbl_ifx_openmpi:0.0) (uses ifx, OpenMPI)

### Obtaining dependencies on supported HPC systems with LMOD

The following module combinations allow ROMS installation on a few tested HPC systems via `module load <package/version>`:

| HPC System    | Modules to Load                                                                                                     | Environment vars                                                                                         |
|---------------|----------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| Derecho       | `ncarenv/23.09`, `intel-oneapi/2024.2.1`, `cray-mpich/8.1.29`, `netcdf/4.9.2`                                        | `MPIHOME=${CRAY_MPICH_PREFIX}/`, `NETCDFHOME=${NETCDF}/`, `LIBRARY_PATH=${LIBRARY_PATH}:{NETCDFHOME}/lib`, `LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDFHOME}/lib` |
| Expanse       | `slurm`, `sdsc`, `DefaultModules`, `shared`, `cpu/0.15.4`, `intel/19.1.1.217`, `mvapich2/2.3.4`, `netcdf-c/4.7.4`, `netcdf-fortran/4.5.3` | `NETCDFHOME=${NETCDF_FORTRANHOME}/`, `MPIHOME=${MVAPICH2HOME}/`, `MPIROOT=${MVAPICH2HOME}/`                                                 |
| Perlmutter    | `cpu/1.0`, `cray-hdf5/1.12.2.9`, `cray-netcdf/4.9.0.9`                                                               | `MPIHOME=${CRAY_MPICH_PREFIX}/`, `NETCDFHOME=${CRAY_NETCDF_PREFIX}/`, `PATH=${PATH}:${NETCDFHOME}/bin`, `LIBRARY_PATH=${LIBRARY_PATH}:${NETCDFHOME}/lib` |


## Obtaining and compiling ROMS
