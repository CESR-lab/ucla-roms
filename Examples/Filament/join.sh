#!/bin/bash

ncjoin --delete L4_ideal_his*.?.nc
ncjoin --delete grid.?.nc
ncview L4_ideal_his.*.nc &

