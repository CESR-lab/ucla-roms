20201120 - DevinD:

do_roms.sh:
This example runs the same Examples/Flux_frc example here with the old lmd_kpp.
It also runs the actual Examples/Flux_frc example with the latest lmd_kpp.
It then does a ncdiff on the results file to see how they diverge.


This is a version of the old code's lmd_kpp, which can give rather different results. 
It has been adapted to work in the new code for comparison.

To use this old version in your roms compilation, add this file lmd_kpp.F to wherever you compile roms,
it will then override the newest lmd_kpp.F.

