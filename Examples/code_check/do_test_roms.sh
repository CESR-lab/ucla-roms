#!/bin/bash

# need this here as well, in case example run on its own:
if [ "$1" != "expanse" -a "$1" != "maya" -a "$1" != "laptop" ]
then
echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
exit
fi

bm_file="benchmark.result_$1"                    # set benchmark specific to machine (maya/expanse)
#echo "$bm_file"

# 1) Run test case:
#echo "  test compiling..."  
cp ../../code_check/diag.F ../                   # need special diagnostic output for full precision
cd ../
make &> /dev/null
echo "  test running..."
mpirun -n 6 ./roms code_check/benchmark.in > code_check/test.log
rm diag.F roms test_his.*.*.nc grid.*.nc &> /dev/null  # grid needed for analytical examples
cd code_check

# 2) Python - confirm values:
python3 test_roms.py $bm_file
retval=$?
#echo $retval

# 3) Rename results logs so they can't be mistakenly read by the 
#    python script even if new simulation doesn't run
mv test.log test_old.log

exit $retval                                     # pass success value onto do_test_all script

# Notes:
# - Compile/ directories are left in the various examples for fast re-compilation of all
#   test cases.






