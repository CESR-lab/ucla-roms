#!/bin/bash
########################################################
#  Define files and run parameters
########################################################
# edit these before using
PROJECTDIR=/u/home/m/minnaho/ucla-roms/Examples/Pipes_ana/
SCRATCHDIR=/u/home/m/minnaho/ucla-roms/Examples/Pipes_ana/
INFILESDIR=${PROJECTDIR}
INFILE=pipes_ana.in

cd $SCRATCHDIR
cp -f $INFILESDIR/$INFILE  .
cp -f $PROJECTDIR/roms .

time mpirun -np 6 -env I_MPI_FABRICS shm:ofa ./roms $INFILE >roms_log.out
