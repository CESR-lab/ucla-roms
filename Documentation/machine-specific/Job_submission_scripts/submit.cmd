#!/bin/csh -f

# replace with your scratch and submit directory here
  set scratchdir = /u/home/m/minnaho/ucla-roms/Examples/Pipes_ana/
  set submitdir =  /u/home/m/minnaho/ucla-roms/Examples/Pipes_ana/

#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  error = Merged with joblog
#$ -o roms.joblog.$JOB_ID
#$ -j y
#  Parallelism:  6-way parallel
#  Resources requested
#$ -pe dc* 6
#$-l h_data=4G,h_rt=1:00:00
#  Name of application for log
#$ -v QQAPP=intelmpi
#  Email address to notify
#$ -M $User@mail # Do not change! Defaults to email on your account
#  Notify at beginning and end of job
# -m bea
# Initialization for mpi parallel execution
#
  unalias *
  set qqversion = 
  set qqapp     = "intelmpi parallel"
  set qqptasks  = 6
  set qqidir    = $submitdir
  set qqjob     = roms
  set qqodir    = $scratchdir
  cd     $submitdir
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "  roms directory:"
  echo "    "$qqidir
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    $qqscratch"
  echo "  roms 6-way parallel job configuration:"
  echo "    $qqconfig" | tr "\\" "\n"
#
  echo ""
  echo "roms started on:   "` hostname -s `
  echo "roms started at:   "` date `
  echo ""
#
# Run the user program
# Load modules

  source /u/local/Modules/default/init/modules.csh
  module purge
  module load intel/2020.4  mpich/3.4
  module load curl/7.70.0
  module load netcdf/c-4.7.4
  module load netcdf/fortran-4.5.3
  module load ncview/2.1.7

  setenv OMP_NUM_THREADS 1

./do_roms_hoffman2.bash

  echo ""
  echo "roms finished at:  "` date `
