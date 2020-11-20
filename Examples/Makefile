# Just type: make 
# (In an example directory or work directory)

# This makefile does:
# 1) Makes a sync of the distribution code to a Compile 
#    directory here (will create directory if not here yet)
# 2) Copies local files in this directory into that Compile directory
# 3) Compiles code into executable for the example. 
# 4) Copies 'roms' from Compile dir into this dir.

# This line de-activates all implicit rules
.SUFFIXES:

.PHONY: all clean depend

all:
	rsync -a $(ROMS_ROOT)/src/*.F $(ROMS_ROOT)/src/*.h  Compile
	@rsync -a $(ROMS_ROOT)/src/Makedefs.inc Compile
	@rsync -a $(ROMS_ROOT)/src/Make.depend Compile
	@rsync -a $(ROMS_ROOT)/src/Makefile Compile
	cp -p *.h *.F Compile 2>/dev/null || :
	cd Compile; make depend; make mpc; make -j6; mv roms ..

# note: "2>/dev/null || :" above is needed to suppress error 
#       stop if no .h or .F files in current directory.

clean:
	rm roms; cd Compile; make clean

# The no-print suppress makes messages of entering/leaving directories
MAKEFLAGS += --no-print-directory
	
