# ifdef LIMIT_UNSTABLE_ONLY
            if (Bfsfc<0.) zscale=min(zscale, hbl(i,j)*epssfc)
# else
            zscale=min(zscale, hbl(i,j)*epssfc)
# endif
# ifdef MASKING
            zscale=zscale*rmask(i,j)
# endif
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3

! Stable regime.

            if (zetahat >= 0.) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5.*zetahat,
     &                                                   1.E-20 )
              ws=wm
            else

! Unstable regime: note that "zetahat" is always negative here,
! also negative are constants "zeta_m" and "zeta_s".

              if (zetahat > zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16.*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat > zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16.*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif

