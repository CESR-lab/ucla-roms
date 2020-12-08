#!/bin/bash

ncjoin --delete pipes_ana_his*.?.nc
ncjoin --delete grid.?.nc
ncview pipes_ana_his.*.nc &

