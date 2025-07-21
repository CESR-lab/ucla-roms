# Obtaining dependencies


ROMS runs on Unix-like systems. It additionally requires:

- A Fortran compiler (`gfortran` and `ifx` are actively supported)
- The `netcdf-fortran` library for netCDF i/o
- An MPI library (e.g. MPICH, OpenMPI, Intel MPI) for parallelization

... Using ROMS with MARBL biogeochemistry additionally requires MARBL v0.45 (see below section on installing MARBL)

This page details different methods for obtaining these dependencies:

```{contents}
:depth: 2
:local:
```

## Using containers
:::{warning}
This is best suited to cloud instances, personal computers, and small-scale scientific computing systems.
On high-performance computing systems, running MPI programs inside containers may be nontrivial. You may wish to consult your HPC's documentation or contact its support team.
:::

The easiest way to obtain a ROMS-ready environment is in a pre-built container. 
We maintain two container images with pre-configured environments containing all of the above dependencies:

- [Image 1](https://ghcr.io/dafyddstephenson/roms_gfortran_build_env:1.0) (uses gfortran, MPICH)
- [Image 2](https://ghcr.io/cworthy-ocean/marbl_ifx_openmpi:0.0) (uses ifx, OpenMPI)

You can fetch and run these images using software like `podman` or `docker`, e.g.:

```
podman pull ghcr.io/cworthy-ocean/marbl_ifx_openmpi:0.0 
podman run -it ghcr.io/cworthy-ocean/marbl_ifx_openmpi:0.0 
```

You can then clone ROMS into the container as described on the {doc}`following page <compilation>`.

## Using conda
:::{warning}
This is best suited to cloud instances, personal computers, and small-scale scientific computing systems.
High-performance computing systems typically have internal package management systems that vary between machines and over time. MPI and other libraries installed via conda may have limited compatibility.
:::

The necessary (GNU) dependencies for ROMS can be installed via conda in an environment created from an included `yaml` file:

```
conda env create -f environments/conda_environment.yml --name roms_env
conda activate roms_env
```

**NOTE** If planning to run ROMS with MARBL biogeochemistry, MARBL will need to be built separately. See the section below on installing MARBL.

You will also need to set some necessary environment variables telling ROMS where to find your conda-installed libraries:

```
export ROMS_ROOT=⚠️ EDIT ME!<your ROMS clone location>⚠️
export MPIHOME=${CONDA_PREFIX}
export NETCDFHOME=${CONDA_PREFIX}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NETCDFHOME/lib" 
export PATH="$ROMS_ROOT/Tools-Roms:$PATH"
```

Verify you set `$ROMS_ROOT` correctly above with `ls $ROMS_ROOT`. You should see directories like `src`,`docs`, etc.
These variables will need to be set every time you wish to use ROMS, so it may be convenient to set their values in a location such as a shell startup file (e.g. `~/.bashrc`)

## Using LMOD on supported HPC systems.

The ROMS repo includes a selection of files for configuring the environment for ROMS on supported HPC systems, in the `environments` directory.
If you are on a supported machine, simply run `source <machinename>.sh` inside that directory (note that paths in these files are relative, so they should not be moved).
If you are running ROMS on another system and would like to share your configuration, or to report that one of these configurations is out of date, please open an [issue](https://github.com/CWorthy-ocean/ucla-roms/issues/new) or [pull request](https://github.com/CWorthy-ocean/ucla-roms/compare).

**NOTE** If planning to run ROMS with MARBL biogeochemistry, MARBL will need to be built separately. See the section below on installing MARBL.

## Installing MARBL

MARBL is an optional dependency used for advanced biogeochemical modeling in ROMS. ROMS works with MARBL v0.45. To install:

1. Clone MARBL to its own directory and checkout the supported version:

```
git clone https://github.com/marbl-ecosys/MARBL.git
cd MARBL/
git checkout marbl0.45.0
```

2. In this directory (the repo top level), set the `MARBL_ROOT` environment variable, so ROMS can find MARBL:

```
export $MARBL_ROOT=$(pwd)
```

This variable will need to be set every time you wish to use MARBL, so it may be convenient to set its value in a location such as a shell startup file (e.g. `~/.bashrc`):


3. Compile MARBL using the same compiler and environment you set up for ROMS above, e.g.

```
cd $MARBL_ROOT/src/
make <compiler> USEMPI=TRUE 
```

4. Where `<compiler>` is either `intel` (if your compiler is `ifx` or `ifort`) or `gnu` (if your compiler is `gfortran`).
If you are unsure, your compiler information should be returned by `mpifort --version`


