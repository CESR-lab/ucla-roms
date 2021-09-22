20201120 - DevinD:

In order to test the merits of the 'first touch' priniciple, we commented out everything in init_arrays.F, except for rmask as the code does not work without it set. We didn't notice any performance change by essentially skipping the first touch, but we need to test it further on more hardware.
