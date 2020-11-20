/* This is include file "bblm.h"
---------------------------------------------------------------------
 Copyright (c) 1996 Rutgers University
---------------------------------------------------------------------
 Ab          Wave bottom excursion amplitude (m).
 Awave       Wind induced wave amplitude (m) at RHO-points.
 Cr          Nondimentional function that determines the importance
                  of currents and wind induced waves on bottom stress
                  at RHO-points.
 Dwave       Wind induced wave direction (radians) at RHO-points.
 Pwave       Wind induced wave period (s) at RHO-points.
 Sdens       Sediment grain density (kg/m^3) at RHO-points.
 Ssize       Sediment grain diameter size (m) at RHO-points.
 Ub          Wave maximum bottom horizontal velocity (m/s).
 UstarC      Time-averaged near-bottom friction current magnitude
                  (m/s) at RHO-points.
----------------------------------------------------------------- */
      real Ab(0:L,0:M),    Awave(0:L,0:M), Cr(0:L,0:M),
     &     Dwave(0:L,0:M), Pwave(0:L,0:M), Sdens(0:L,0:M),
     &     Ssize(0:L,0:M), Ub(0:L,0:M),    UstarC(0:L,0:M)
      common /fbblm/ Ab, Awave, Cr, Dwave, Pwave,
     &                   Sdens, Ssize, Ub, UstarC
 
