#!/bin/bash

ncjoin --delete ideal_his*.?.nc
ncjoin --delete grid.?.nc
ncview ideal_his.*.nc &

