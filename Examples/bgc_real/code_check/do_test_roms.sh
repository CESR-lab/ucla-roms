#!/bin/bash

# need this here as well, in case example run on its own:
if [ "$1" != "expanse" -a "$1" != "maya" -a "$1" != "laptop" -a "$1" != "github" ]
then
echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
exit
fi

bm_file="benchmark.result_$1"                    # set benchmark specific to machine (maya/expanse)
echo "$bm_file"
retval=0
for BGC_MODEL in {"BEC","MARBL"};do
    echo "Running bgc_real test with ${BGC_MODEL}"
    # 1) Compile test case:
    echo "  test compiling..."    

    cp -p ../*.h . &> /dev/null
    cp -p ../cppdefs_${BGC_MODEL}.opt ./cppdefs.opt
    cp -p ../param.opt .
    if [ ${BGC_MODEL} == "MARBL" ];then
	cp -p ../marbl_in .
    fi
    cp -p $ROMS_ROOT/Examples/code_check/diag.opt .
    cp -p $ROMS_ROOT/Examples/code_check/Makedefs.inc .
    cp -p $ROMS_ROOT/Examples/Makefile .
    cp -p $ROMS_ROOT/Examples/code_check/test_roms.py .
    make compile_clean &> /dev/null
    make &> /dev/null


    # 2) Run test case:
    echo "  test running..."
    if [ "$1" = "expanse" ]
    then
	srun --mpi=pmi2 -n 6 ./roms benchmark.in > test.log
    else
	mpirun -n 6 ./roms benchmark.in > test.log
    fi

    rm *.h       &> /dev/null
    rm *.nc      &> /dev/null
    rm diag.opt  &> /dev/null
    rm Make*     &> /dev/null
    rm param.opt &> /dev/null
    rm cppdefs.opt &> /dev/null
    rm roms      &> /dev/null
    if [ ${BGC_MODEL} == "MARBL" ];then
	rm marbl_*   &> /dev/null
    fi
    
    # 2) Python - confirm values:
    python3 test_roms.py $bm_file
    retval_tmp=$?
    echo "exit code for ${BGC_MODEL} test: $retval_tmp"
    if [ $retval_tmp -gt $retval ];then
	echo "Test failed for ${BGC_MODEL}"
	retval=retval_tmp
    fi

    # 3) Rename results logs so they can't be mistakenly read by the 
    #    python script even if new simulation doesn't run
    mv test.log test_old_${BGC_MODEL}.log

    rm test_roms.py
done
exit $retval                                     # pass success value onto do_test_all script

# Notes:
# Compile directories are left in the various examples for fast re-compilation of all
#  test cases.






