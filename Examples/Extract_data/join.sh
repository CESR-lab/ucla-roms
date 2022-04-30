#!/bin/bash

ncjoin -d sample_his*.?.nc
extract_data_join sample_ext.00000.?.nc
ncview grid1_bry.00000.nc &
