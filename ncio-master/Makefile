
FC = ifort
DFLAGS = -O3


## Shared library 
lib/libncio.so: ncio.f90
	$(FC) -O3 -c ncio.f90 -module include
	xiar rc lib/libncio.a ncio.o
#$(FC) -O3 -c -shared -fPIC -o lib/libncio.so ncio.f90 -module include


clean:
	rm lib/* include/* *.o *.mod

