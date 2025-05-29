#!/bin/bash

# need this here as well, in case example run on its own:
case "$1" in
    expanse|maya|laptop|github|github_ifx)
	echo "running test for $1"
	;;
    *)
	echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
	exit
	;;
esac
    
bm_file="benchmark.result_$1"                    # set benchmark specific to machine (maya/expanse)
echo "$bm_file"

# 1) Compile test case:
echo "##############################"
echo "  test compiling..."
echo "##############################"

cp -p ../*.h . &> /dev/null
cp -p ../cppdefs.opt .
cp -p ../param.opt .
cp -p $ROMS_ROOT/Examples/code_check/diag.opt .
#cp -p $ROMS_ROOT/Examples/code_check/Makedefs.inc .
cp -p $ROMS_ROOT/Examples/Makefile .
make compile_clean &> /dev/null
make > compile.log 


# 2) Run test case:
echo "##############################"
echo "  test running..."
echo "##############################"

if [ "$1" = "expanse" ]
then
    srun --mpi=pmi2 -n 6 ./roms benchmark.in > test.log
else
    mpirun -n 6 ./roms benchmark.in > test.log #2>&1 | tee -a test.log
fi


rm *.h       &> /dev/null
rm *.nc      &> /dev/null
rm diag.opt  &> /dev/null
rm Make*     &> /dev/null
rm param.opt &> /dev/null
rm cppdefs.opt# &> /dev/null
rm roms      &> /dev/null

# 2) Python - confirm values:
cp $ROMS_ROOT/Examples/code_check/test_roms.py . 
python3 test_roms.py $bm_file
retval=$?

rm test_roms.py

# 3) Rename results logs so they can't be mistakenly read by the 
#    python script even if new simulation doesn't run
mv test.log test_old.log

exit $retval                                     # pass success value onto do_test_all script

# Notes:
# Compile directories are left in the various examples for fast re-compilation of all
#  test cases.






