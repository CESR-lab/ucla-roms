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
	rsync -r $(ROMS_ROOT)/Makedefs Compile
	rsync -r $(ROMS_ROOT)/netcdf.inc Compile
	rsync    $(ROMS_ROOT)/Make.depend Compile
	rsync    $(ROMS_ROOT)/Makefile Compile
	cp -p *h *F Compile
	cd Compile; make -j8; mv roms ../roms_tank


	
