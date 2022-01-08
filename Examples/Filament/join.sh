#!/bin/bash

rm grid.*
ncjoin --delete ideal_his*.?.nc
ncview ideal_his.*.nc &

