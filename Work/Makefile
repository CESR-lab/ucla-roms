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

optfiles := $(shell ls *.opt)
tag := $(shell git --git-dir $(ROMS_ROOT)/.git rev-parse HEAD)
tline1 := git_hash=
tline2:= 'git_hash=\"$(tag)\"'

all:
	rsync -a $(ROMS_ROOT)/src/*.F $(ROMS_ROOT)/src/*.h  Compile
	make show_tag
	@rsync -a $(ROMS_ROOT)/src/*.opt Compile	
#@rsync -a $(ROMS_ROOT)/src/Makedefs.inc Compile
	@rsync -a $(ROMS_ROOT)/src/Make.depend Compile
	@rsync -a $(ROMS_ROOT)/src/Makefile Compile
	cp -p *.h *.F *.opt Compile 2>/dev/null || :
	cd Compile; make depend 2>/dev/null || :
	cd Compile; make -j6; mv roms ..

# note: "2>/dev/null || :" above is needed to suppress error 
#       stop if no .h or .F files in current directory.
#       Also, makedepf90 in make depend can throw out strange
#       errors for .h files for inexplicable changes made to
#       cppdefs.h.


show_tag:
	@sed -i '/$(tline1)/c\      $(tline2)' Compile/add_git_hash.F
compile_clean:
	rm -r Compile/ roms 2>/dev/null || :
	
work_clean:
	rm -r Compile/ *.F *.h *.in *.sh *.nc roms 2>/dev/null || :	

code_check_clean:
	@cd $(ROMS_ROOT)/Examples/WEC_real   ; make compile_clean
	@cd $(ROMS_ROOT)/Examples/Rivers_ana ; make compile_clean
	
copy_to_Work_dir:
	cp -p *.F *.h *.in *.sh ../../Work 2>/dev/null || :

nhmg:
	cd $(ROMS_ROOT)/NHMG/src; make clean; cd .. ; make		

# The no-print suppress makes messages of entering/leaving directories
MAKEFLAGS += --no-print-directory
	
