# This makefile make a sync of the distribution code to a compile directory,
# After that, it copies local file into that directory
# and make the executable for the example.
# Just type: make

# This line de-activates all implicit rules
.SUFFIXES:

.PHONY: all clean depend

all:
#rsync -a $(ROMS_ROOT)/*.F $(ROMS_ROOT)/*.h $(ROMS_ROOT)/*.mod Compile
	rsync -a $(ROMS_ROOT)/*.F $(ROMS_ROOT)/*.h  Compile
	rsync -a $(ROMS_ROOT)/Makedefs.inc Compile
	rsync -a $(ROMS_ROOT)/netcdf.inc Compile
	rsync -a $(ROMS_ROOT)/Make.depend Compile
	rsync -a $(ROMS_ROOT)/Makefile Compile
	cp -p *.h  Compile
	cd Compile; make mpc; make -j8; mv roms ..


	
