#!/bin/bash

# need this here as well, in case example run on its own:
case "$1" in
    expanse|maya|laptop|github_gnu|github_ifx)
	echo "running test for $1"
	;;
    *)
	echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
	exit 1
	;;
esac

case "$2" in
    MARBL|BEC)
	echo "running test for $2"
	;;
    *)
	echo "Script must have second argument 'MARBL' or 'BEC'."
	;;
esac
BGC_MODEL=$2

bm_file="benchmark.result_$1"                    # set benchmark specific to machine (maya/expanse)
echo "$bm_file"

echo "Running bgc_real test with ${BGC_MODEL}"
# 1) Compile test case:
echo "##############################"
echo "  test compiling, ${BGC_MODEL}..."
echo "##############################"

cp -p ../*.h . &> /dev/null
cp -p ../cppdefs_${BGC_MODEL}.opt ./cppdefs.opt
cp -p ../param.opt .
if [ ${BGC_MODEL} == "MARBL" ];then
    cp -p ../marbl_in .
fi
cp -p $ROMS_ROOT/Examples/code_check/diag.opt .
#cp -p $ROMS_ROOT/Examples/code_check/Makedefs.inc .
cp -p $ROMS_ROOT/Examples/Makefile .
cp -p $ROMS_ROOT/Examples/code_check/test_roms.py .
make compile_clean &> /dev/null
make > compile_${BGC_MODEL}.log # 2>&1 | tee -a compile_${BGC_MODEL}.log

# 2) Run test case:
echo "##############################"
echo "  test running, ${BGC_MODEL} ..."
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
rm cppdefs.opt &> /dev/null
rm roms      &> /dev/null

# 2) Python - confirm values:
python3 test_roms.py $bm_file
retval=$?
echo "exit code for ${BGC_MODEL} test: $retval"

# 3) Rename results logs so they can't be mistakenly read by the
#    python script even if new simulation doesn't run
mv test.log test_old.log
cp test_old.log test_old_${BGC_MODEL}.log

rm test_roms.py
exit $retval                                     # pass success value onto do_test_all script

# Notes:
# Compile directories are left in the various examples for fast re-compilation of all
#  test cases.






