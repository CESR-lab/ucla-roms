#!/bin/bash

ncjoin -d sample_his*.?.nc
ncjoin -d sample_dia.*.nc
ncview sample_dia.*.nc
