#!/bin/bash

ncjoin --delete pipes_grd*.?.nc
ncjoin --delete pipes_his*.?.nc
ncview pipes_his.*.nc &

