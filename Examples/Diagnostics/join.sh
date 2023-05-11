#!/bin/bash

ncjoin -d diags_his*.?.nc
ncjoin -d diags_dia*.?.nc
ncview diags_dia.*.nc
