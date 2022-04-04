#!/bin/bash

ncjoin -d sample_his*.?.nc
join_child_bry -d sample_child_bry.00000.?.nc
join_child_bry -d sample_child_bry_2.00000.?.nc
join_child_bry -d sample_child_bry.00002.?.nc
join_child_bry -d sample_child_bry_2.00002.?.nc
ncview sample_child_bry.0000.nc &
