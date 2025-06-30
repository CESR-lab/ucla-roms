#!/bin/bash

# ******** USER INPUT START ************
# declare an array of Example folders and .in names to use:
declare -a Examples=( "Flux_frc" "Pipes_ana" "Pipes_real" "Rivers_ana" "Rivers_real" "Filament" "bgc_real" )
# ******** USER INPUT END   ************

arg=$1

case "$arg" in
    expanse|maya|laptop|github_gnu|github_ifx)
	echo "running test for $arg"
	;;
    *)
	echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
	exit 1
	;;
esac

error_cnt=0                                      # count of exit codes from each test
total=${#Examples[*]}                            # total number of examples
for (( i=0; i<=$(( $total -1 )); i++ ))          # run test cases:
do
    cd ../${Examples[i]}/code_check/
    echo "##############################"
    echo "${Examples[i]}:"
    echo "##############################"

    case ${Examples[i]} in
	bgc_real)
	    ./do_test_roms.sh $arg BEC
	    retval=$?
	    error_cnt=$(( $error_cnt + $retval ))
	    ./do_test_roms.sh $arg MARBL
	    ;;
	*)
	    ./do_test_roms.sh $arg
	    ;;
    esac
    retval=$?                                      # $? gives exit code from ./do_test_roms.sh
    error_cnt=$(( $error_cnt + $retval ))

  if [ $retval -ne 0 ]
  then
    echo -e "  test failed! \n"
#   break
  fi

  cd ../                                         # need return out of /code_check/ for next iteration
done

if [ $error_cnt -eq 0 ]
then
    echo "ALL TESTS SUCCESSFUL!"
    exit 0
else
    echo "ERROR - A TEST FAILED!"
    exit 1
fi

# Notes:
# - Compile/ directories are left in the various examples for fast re-compilation of all
#   test cases.






