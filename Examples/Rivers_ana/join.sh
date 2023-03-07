#!/bin/bash

ncjoin --delete river_his*.?.nc
ncjoin --delete river_grd.?.nc
ncview river_his.*.nc

