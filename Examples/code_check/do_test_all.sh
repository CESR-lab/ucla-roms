#!/bin/bash

# ******** USER INPUT START ************
# declare an array of Example folders and .in names to use:
#declare -a Examples=( "Flux_frc" "Pipes_ana" "Pipes_real" "Rivers_ana" "Rivers_real" "Tracers_passive" "WEC_real" )
declare -a Examples=( "Flux_frc" "Pipes_ana" "Pipes_real" "Rivers_ana" "Rivers_real" "Filament" "bgc_real" )
# ******** USER INPUT END   ************

arg=$1

if [ "$arg" != "expanse" -a "$arg" != "maya" -a "$arg" != "laptop" -a "$arg" != "github" ]
then
echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
exit
fi
                                                 
error_cnt=0                                      # count of exit codes from each test
total=${#Examples[*]}                            # total number of examples
for (( i=0; i<=$(( $total -1 )); i++ ))          # run test cases:
do
  cd ../${Examples[i]}/code_check/
  echo "${Examples[i]}:"
  ./do_test_roms.sh $arg
  
  retval=$?                                      # $? gives exit code from ./do_test_roms.sh
  error_cnt=$(( $error_cnt + $retval ))          
  if [ $error_cnt -gt 0 ]
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






