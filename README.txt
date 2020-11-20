ROMS basic usage:
-----------------

1) Copy files you need to physically edit from src/ into:
   A) Work/ if you are working on your own simulation config,
   or
   B) Example/<example>/ if you want to run a pre-existing example.
   
2) Edit your files in A) or B) as applicable.

3) Hit 'make' in directory A) or B) to compile roms.

   For more info see readme's in Documentation/ directory.
   Do not edit code in the src/ directory!

Directories:
------------

Relevant:

   Documentation/ - Your next step for more info...

   src/           - The latest roms code is kept here.
                    Do not modify code here! use Work/ or Examples/...

   Examples/      - Collection of previously run examples, however
                    without netcdf files. Stored elsewhere.
                    
   Tools-Roms/    - Useful pre & post processing tools can be compiled here.
   
   Work/          - Edit roms files here and compile new simulations here.                
                   
Less relevant: (you're unlikely to use these)                   
                 
   Legacy/        - Old code no longer used in latest version of ROMS.
      
   Makedefs/      - Machine files for running on different machines
   
   NHMG/          - Non-hydrostatic modules   
                                   



Examples/ - 
