#!/bin/bash

# ******** USER INPUT START ************
# Declare an array of Example folders and .in names to use:
#                      TEST 1         TEST 2        TEST 3
declare -a Examples=("WEC_real"     "Rivers_ana"  "WEC_real"    )
#   results folder in benchmarks/
declare -a      dir=("wec"          "river"       "netcdf"      )
#   name of .in file without .in extension
declare -a     file=("sample_wec"   "river_ana"   "wec_netcdf"  )
# ******** USER INPUT END **************

# 1) Run test cases: first copy files to this directory, run simulation, delete files leaving netcdf.
total=${#Examples[*]} # Total number of examples
for (( i=0; i<=$(( $total -1 )); i++ ))
do
  cd ../${Examples[i]}
  echo "${Examples[i]} test compiling..."
  make &> /dev/null
  mv roms ../code_check
  cd ../code_check
  # .in file needs to be in this directory as can't have more than 3 leading periods ../../../ not allowed.
  cp benchmarks/${dir[i]}/${file[i]}_benchmark.in .
  ./benchmarks/${dir[i]}/do_BM.sh # run roms - output piped to benchmarks/wec/*.log
  rm ${file[i]}_benchmark.in roms
done

# 2) Python - confirm values:
python3 test_roms.py

# 3) Rename results logs so they can't be mistakenly read by python 
# script even if new simulation doesn't run
for (( i=0; i<=$(( $total -1 )); i++ ))
do
  mv benchmarks/${dir[i]}/${file[i]}_test.log benchmarks/${dir[i]}/${file[i]}_test_old.log
done

# Notes:
# - Compile/ directories are left in the various examples for fast re-compilation of all
#   test cases.






