#!/bin/bash

# ******** USER INPUT START ************
# declare an array of Example folders and .in names to use:
#                      TEST 1      TEST 2      TEST 3
declare -a Examples=( "Flux_frc" "Pipes_ana" "Pipes_real" "Rivers_ana" "Rivers_real" "Tracers_passive" "WEC_real" )
# ******** USER INPUT END   ************
                                                 
error_cnt=0                                      # count of exit codes from each test
total=${#Examples[*]}                            # total number of examples
for (( i=0; i<=$(( $total -1 )); i++ ))          # run test cases:
do
  cd ../${Examples[i]}/code_check/
  echo "${Examples[i]} test compiling..."
  ./do_test_roms.sh
  
  retval=$?                                      # $? gives exit code from ./do_test_roms.sh
  error_cnt=$(( $error_cnt + $retval ))          
  if [ $error_cnt -gt 0 ]
  then
    echo "  test failed!"
    break
  fi  
  
  cd ../                                         # need return out of /code_check/ for next iteration
done

if [ $error_cnt -eq 0 ]
then
  echo "ALL TESTS SUCCESSFUL!"
else
  echo "ABORTED - A TEST FAILED!"
fi

# Notes:
# - Compile/ directories are left in the various examples for fast re-compilation of all
#   test cases.






