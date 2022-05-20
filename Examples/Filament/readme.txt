2022/01:
Best results achieved with NO initial along filament flow (u0=0).
Therefore online analytical version is preferred over python generated initial conditions,
as there is no need for offline TTW sovler.
LIN_EOS essential to avoid large start-up waves. This is because analytical temperature gives
inital rho, only easily calculated with LIN_EOS.

2021/11:
Idealised_Filament python scripts created by DanielD. It includes TTW solver.
Run script in Idealized_Filament/
python3 setup_grid_init_bry.py params_fil.py

Edit shape of filament by changing parameters in params_fil.py. Grid can be changed there to give desired
grid size for ROMS simulation.

Scritps will produce grid file and initial file for simulation.

ROMS analytical scripts here create geostrophic flow only for comparison.
