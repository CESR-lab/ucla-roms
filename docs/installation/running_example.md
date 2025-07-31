# Compiling and running an example configuration

## Compiling
ROMS does not have a single, central executable, and each configuration must be compiled individually. This is because ROMS takes advantage of several configuration-dependent compile-time optimizations. ROMS includes a number of simple example configurations in the `Examples` directory that can easily be compiled and run on a modern laptop. To verify everything has gone as expected up to this point, we'll compile the `Rivers_real` example, which is a small simulation of a river outflow using netCDF input and output on 6 processors.

First, obtain the necessary netCDF input files for the examples:

```
cd $ROMS_ROOT/Examples/input_data/
./get_input_files.sh
```

Next, compile the model:

```
cd $ROMS_ROOT/Examples/Rivers_real
make
```

You should now have an executable, `roms`, in the directory.

## Running
To run the model, do

```
mpirun -n 6 ./roms rivers.in
```

where `-n 6` is the number of processors, and `rivers.in` is a text file defining runtime parameters such as the number of time steps.
Unless edited, the model will run for 50 time steps.

:::{warning}
On HPC systems, you will likely need to wrap the above MPI run command in a script and submit it to a scheduler. Consult your system's documentation.
:::

The screen should fill with live output from the model.


## Processing the output
Once the run has finished, the directory will contain six files like `rivers_his.20121209133435.0.nc`, which are the model output. The six files each correspond to a different processor, each handling one sixth of the domain. We need to stitch these files together to view the domain in its entirety. This can be done using the `ncjoin` tool included with ROMS:

```
ncjoin rivers_his.20121209133435.*.nc
```

after which there will be an additional netCDF file `rivers_his.20121209133435.nc` containing the joined outputs from the 6 processors - this output file is now ready for analysis.
