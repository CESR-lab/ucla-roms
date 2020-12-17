bgc: porting info from ETH bec2 functionality into roms code. 20201216.

Name changes only for clarity/grouping of bgc code:
ecosys_bec2.h -> bgc_ecosys_bec2.h
bio_diag.h -> bgc_bio_diag.h
tracers.h ->  bgc_tracers.h (though indxPO4, etc is not used anymore by new code for tracer netcdf output, it has been included in order to get ETH port to compile only)
def_his.F -> def_his_bgc_diag() in bgc.F
          -> def_avg_bgc_diag() in bgc.F
def_rst.F -> def_rst_bgc_diag() in bgc.F    
wrt_his.F -> wrt_his_bgc_diag() in bgc.F
wrt_avg.F -> wrt_avg_bgc_diag() in bgc.F
set_avg.F -> set_avg_bgc_diag() in bgc.F
ecosys_bec2.F -> bgc_ecosys_bec2.F
ecosys_bec2_init.F -> bgc_ecosys_bec2_init.F

init_scalars_bec2.F -> init_scalars_bec2() in bgc.F      

Code moved to different files:
iPO4, iNH4, etc in param.h now moved to tracers.F module through include files tracers_defs_idx.h and tracers_defs.h.

Issues:
neither ETH code nor new code works with just Ncycle_SY flag for ncycle.
