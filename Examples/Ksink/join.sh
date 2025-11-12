#!/bin/bash

ncjoin --delete ksink_his*.?.nc
ncjoin --delete ksink_dia*.?.nc
ncview ksink_his.*.nc &
