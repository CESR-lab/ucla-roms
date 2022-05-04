#!/bin/bash

# ******** USER INPUT START ************
# declare an array of Example folders and .in names to use:
declare -a Examples=( "Flux_frc" "Pipes_ana" "Pipes_real" "Rivers_ana" "Rivers_real" "Tracers_passive" "WEC_real" )
# ******** USER INPUT END   ************

arg=$1

if [ "$arg" != "expanse" -a "$arg" != "maya" -a "$arg" != "laptop" ]
then
echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
exit
fi

bm_file="benchmark.result_$arg"                    # set benchmark specific to machine (maya/expanse)

total=${#Examples[*]}                            # total number of examples
for (( i=0; i<=$(( $total -1 )); i++ ))          # run test cases:
do
  cd ../${Examples[i]}/code_check/
  echo "${Examples[i]} updating benchmark..."
  cp test_old.log $bm_file
    
  
  cd ../                                         # need return out of /code_check/ for next iteration
done


echo -e "\nBENCHMARK UPDATE COMPLETE!"






