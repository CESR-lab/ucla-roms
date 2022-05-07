#!/bin/bash

ncjoin --delete pipes_ana_his*.?.nc
rm grid.*.nc
ncview pipes_ana_his.*.nc &

