/* Buffer array to allow reshaping input/output data for compatibility
 with the plotting package and other pre- and post-processing software.
 This implies that variables defined at RHO-, U-, V-, and PSI-points
 when written into netCDF files are dimensioned as follows:
 
 Location     name                  dimensions
 
   RHO-        r2dvar   zeta-type   (0:Lm+1,0:Mm+1)
   PSI-        p2dvar   vort-type   (1:Lm+1,1:Mm+1)
   U-          u2dvar   ubar-type   (1:Lm+1,0:Mm+1)
   V-          v2dvar   vbar-type   (0:Lm+1,1:Mm+1)
 
   RHO-,RHO-   r3dvar   RHO-type    (0:Lm+1,0:Mm+1,  N)
   PSI-,PSI-   p3dvar   PSI-type    (1:Lm+1,1:Mm+1,  N)
   U-,RHO-     u3dvar   U-type      (1:Lm+1,0:Mm+1,  N)
   V-,RHO-     v3dvar   V-type      (0:Lm+1,1:Mm+1,  N)
   RHO-,W-     w3dvar   W-type      (0:Lm+1,0:Mm+1,0:N)
*/
      real buff((Lm+5)*(Mm+5)*(N+1))
      common /buffer/ buff
 
 
