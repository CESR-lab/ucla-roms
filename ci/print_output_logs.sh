#!/bin/bash

# 

# Get the "Examples" array from Examples/code_check/do_test_all.sh
{ 
  IFS= read -r shebang_line 
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" == *"USER INPUT END"* ]] && break
    eval "$line"
  done 
} < do_test_all.sh

if [ ${#Examples[@]} -eq 0 ]; then
    echo "The Examples array is empty."
    echo "This script expects a bash array `Examples` to be defined in"
    echo "Examples/code_check/do_test_all.sh"
    echo "above the line containing `USER INPUT END`."
    echo "Please verify the contents of this file."
    exit 1
fi

# Loop over the "Examples" array
for example in "${Examples[@]}"; do
    echo "----------------------------------------"
    echo "$example"
    echo "----------------------------------------"
    cat ${ROMS_ROOT}/Examples/${example}/code_check/test_old.log
done
