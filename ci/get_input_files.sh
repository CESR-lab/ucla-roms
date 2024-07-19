#!/bin/bash

if command -v curl &> /dev/null;then
    DOWNLOAD_COMMAND="curl -O -L"
elif command -v wget &> /dev/null;then
    DOWNLOAD_COMMAND=wget
fi
echo ${DOWNLOAD_COMMAND}   
# URLs to download
URL_PREFIX="https://github.com/dafyddstephenson/ucla_roms_examples_input_data/raw/main"
files=(
    "roms_LFfrc.nc"
    "roms_bgcflx.nc"
    "roms_default_bgcbry.nc"
    "roms_init_trace.nc"
    "sample_grd_riv.nc"
    "sample_rad_units_DPD.nc"
    "sample_tra_units_DPD.nc"
    "sample_wnd.nc"
    "roms_bgcbry.nc"
    "roms_bry_trace.nc"
    "roms_init_bgc.nc"
    "sample_flux_frc.nc"
    "sample_prec_units_DPD.nc"
    "sample_tides.nc"
    "sample_tracers.nc"
    "sample_wwv_riv.nc"
)

for fname in "${files[@]}";do
    echo "#######################################################"
    echo "FETCHING FILE ${fname}"
    echo "#######################################################"
    ${DOWNLOAD_COMMAND} "${URL_PREFIX}/${fname}"
    partit 3 2 "${fname}"
    #echo "${URL_PREFIX}/${fname}"
done

