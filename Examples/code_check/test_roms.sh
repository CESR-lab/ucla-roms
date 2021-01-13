#!/bin/bash

# 1) Run test cases: first copy files to this directory, run simulation, delete files leaving netcdf.

# SET UP AS LOOP TO ADD MORE EXAMPLES!!

# - Declare an array of Example folders and .in names to use:
#declare -a folders=("WEC_real" "Rivers_ana")

# - a) WEC test:
cd ../WEC_real
#cd ../$folders[0]
echo 'WEC_real test compiling...'
make &> /dev/null
mv roms ../code_check
cd ../code_check
  # .in file needs to be in this directory as can't have more than 3 leading periods ../../../ not allowed.
cp benchmarks/wec/sample_wec_benchmark.in .
./benchmarks/wec/do_BM.sh # run roms - output piped to benchmarks/wec/*.log
rm sample_wec_benchmark.in roms

# - b) Another test:
cd ../Rivers_ana
echo 'Rivers_ana test compiling...'
make &> /dev/null
mv roms ../code_check
cd ../code_check
  # .in file needs to be in this directory as can't have more than 3 leading periods ../../../ not allowed.
cp benchmarks/river/river_ana_benchmark.in .
./benchmarks/river/do_BM.sh
rm river_ana_benchmark.in roms grid.?.nc

# 2) Python - confirm values:
python3 test_roms.py

# 3) Rename results logs so they can't be mistakenly read by python 
# script even if new simulation doesn't run
mv benchmarks/wec/sample_wec_test.log benchmarks/wec/sample_wec_test_old.log
mv benchmarks/river/river_ana_test.log benchmarks/river/river_ana_test_old.log

# Notes:
# - Compile/ directories are left in the various examples for fast re-compilation of all
#   test cases.






