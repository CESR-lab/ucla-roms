#!/bin/bash

#mpiexec -np 6 roms river_ana.in < /dev/null > & jobout &
mpiexec -np 1 ./roms roms.in 
