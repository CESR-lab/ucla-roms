#!/bin/bash

ncjoin --delete river_his*.?.nc
ncjoin --delete river_grid.?.nc
ncview river_ana_his.*.nc

