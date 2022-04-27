#!/bin/bash

ncjoin -d sample_his*.?.nc
extract_data_join sample_ext.00000.?.nc
ncview grid1.00000.nc &
