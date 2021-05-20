#!/bin/bash

ncjoin -d sample_his*.?.nc
ncjoin -d sample_on_diags.*.nc
ncview sample_on_diags.*.nc
