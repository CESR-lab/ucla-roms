#!/bin/bash
# Run the checkerboard test with the Dirichlet and the Neumann surface condition
# for the NH pressure (same executable, namelist switch), join the output and
# print the growth table for each.  Usage: make && ./run_both.sh
set -e
for bc in dirichlet neumann; do
  d=run_$bc; rm -rf $d; mkdir $d; cd $d
  cp ../roms ../tank.in .
  if [ $bc = neumann ]; then sed 's/surface_neumann *= *\.false\./surface_neumann = .true./' ../nhmg_namelist > nhmg_namelist
  else cp ../nhmg_namelist .; fi
  mpiexec -np 4 ./roms tank.in > run.log 2>&1
  for base in $(ls tank_*.?.nc | sed -E 's/\.[0-9]\.nc$//' | sort -u); do ncjoin -d $base.?.nc > /dev/null; done
  echo "== $bc"; python3 ../checker.py . | awk 'NR==1 || $1==0.0 || $1==1.0 || $1==2.0 || $1==3.0 || $1==5.0 || $1==10.0 || $1==20.0 || $1==50.0'
  cd ..
done
