      subroutine lmd_kpp_tile (istr,iend,jstr,jend,       Kv,Kt,Ks,
     &                         ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,
     &                         Cr,FC, wrk1,wrk2,
     &                         Gm1,dGm1dS,  Gt1,dGt1dS,  Gs1,dGs1dS,
     &                                                   kbls, kmo)
      implicit none
      integer(kind=4), parameter ::
     &              LLm=1504,  MMm=1024 ,N=50
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , NT=2
      integer(kind=4) istr,iend,jstr,jend, i,j,k
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,0:N) :: Kv, 
     &                               Kt, Ks
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: ustar, 
     &                             Bo, Bosol
     &                                                  , hbl, bbl
     &                                              , FX, FE, FE1
      real(kind=8), dimension(istr-2:iend+2,0:N) :: Cr,FC, wrk1,wrk2
      real(kind=8), dimension(istr-2:iend+2) ::  Bfsfc_bl,
     &                               Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS
      integer(kind=4), dimension(istr-2:iend+2) :: kbls, kmo
      real(kind=8), parameter ::
     &   Ricr=0.45_8,
     &   Ri_inv=1._8/Ricr,
     &   epssfc=0.1_8,
     &   betaT=-0.2_8,
     &   nubl=0.01_8,
     &   nu0c=0.1_8,
     &   Cv=1.8_8,
     &   C_MO=1._8,
     &   C_Ek=258._8,
     &   Cstar=10._8,
     &   zeta_m=-0.2_8,
     &   a_m=1.257_8,
     &   c_m=8.360_8,
     &   zeta_s=-1.0_8,
     &   a_s=-28.86_8,
     &   c_s=98.96_8
      real(kind=8), parameter :: r2=0.5_8, r3=1._8/3._8, r4=0.25_8, 
     &                             EPS=1.D-20
      real(kind=8) Cg, ustar3, Bfsfc, zscale, zetahat, ws,wm, Kern, 
     &                             Vtc,Vtsq,
     &         ssgm, z_bl, Av_bl,dAv_bl, At_bl,dAt_bl,  As_bl,dAs_bl,
     &                         f1,a1,a2,a3,  cff,cff1, cff_up,cff_dn
      real(kind=8) h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_h/h /grd_hinv/hinv /grd_f/f /grd_fomn/fomn
      real(kind=8) angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_angler/angler
      real(kind=8) latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_latr/latr /grd_lonr/lonr
      real(kind=8) pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_pm/pm     /mtrix_pn/pn
     &       /mtrix_dm_r/dm_r /mtrix_dn_r/dn_r
     &       /mtrix_dm_u/dm_u /mtrix_dn_u/dn_u
     &       /mtrix_dm_v/dm_v /mtrix_dn_v/dn_v
     &       /mtrix_dm_p/dm_p /mtrix_dn_p/dn_p
      real(kind=8) iA_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) iA_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_iAu/iA_u  /mtrix_iAv/iA_v
      real(kind=8) dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_dmde/dmde   /mtrix_dndx/dndx
      real(kind=8) pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_pmon_u/pmon_u /mtrix_pnom_v/pnom_v
     &                            /mtrix_grdscl/grdscl
      real(kind=8) rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask /mask_p/pmask
     &       /mask_u/umask /mask_v/vmask
      real(kind=8) u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real(kind=8) v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real(kind=8) t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real(kind=8) FlxU(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) FlxV(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real(kind=8) Wi(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /flx_FU/FlxU /flx_FV/FlxV /flx_We/We /flx_Wi/Wi
      real(kind=8) Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_zw/z_w /grid_zr/z_r /grid_Hz/Hz
      real(kind=8) sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_sustr/sustr /frc_svstr/svstr
      real(kind=8) srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_srflx/srflx
      real(kind=8) stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /frc_stflx/stflx
      real(kind=8) visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      real(kind=8) diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2/diff2
      real(kind=8) Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real(kind=8) Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,isalt)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real(kind=8) bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real(kind=8) hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /kpp_hbl/hbls
      real(kind=8) swr_frac(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /kpp_swr_frac/swr_frac
      real(kind=8) ghat(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /kpp_ghat/ghat
      real(kind=4) cpu_init, cpu_net
      real(kind=8) WallClock, time, tdays
      integer(kind=4) proc(2), numthreads, iic, kstp, knew
     &                           , iif, nstp, nnew, nrhs
     &                           , priv_count(16)
      logical synchro_flag, diag_sync
      common /priv_scalars/ WallClock, cpu_init, cpu_net,
     &   proc, time, tdays, numthreads, iic,  kstp, knew
     &                           , iif, nstp, nnew, nrhs
     &       , priv_count, synchro_flag, diag_sync
C$OMP THREADPRIVATE(/priv_scalars/)
      real(kind=8) start_time, dt, dtfast, time_avg
      real(kind=8) xl,el, rdrg,rdrg2,Zob, visc2,gamma2
      common /scalars_main/ start_time, dt, dtfast, time_avg,
     &                   xl,el, rdrg,rdrg2,Zob, visc2,gamma2
      real(kind=8) rho0, tnu2(NT)
      common /scalars_main/ rho0, tnu2
      real(kind=8) dSSSdt
      common /scalars_main/ dSSSdt
      real(kind=8) v_sponge
      common /scalars_main/ v_sponge
      real(kind=8) ubind
      common /scalars_main/ ubind
      integer(kind=4) ntstart, ntimes, ndtfast, nfast, ninfo, 
     &                           may_day_flag,
     &                                                barr_count(16)
      common /scalars_main/ ntstart, ntimes, ndtfast, nfast, ninfo,
     &                               may_day_flag,    barr_count
      integer(kind=4) forw_start
      common /scalars_main/ forw_start
      real(kind=8), parameter :: pi=3.14159265358979323_8, 
     &                        Eradius=6371315._8,
     &              deg2rad=pi/180._8, rad2deg=180._8/pi, 
     &                         day2sec=86400._8,
     &                   sec2day=1._8/86400._8, Cp=3985._8, 
     &                           vonKar=0.41_8
     &                 , g=9.81_8
      integer(kind=4) imin,imax,jmin,jmax
      if (istr==iwest .and. .not.west_exchng) then
        imin=istr
      else
        imin=istr-1
      endif
      if (iend==ieast .and. .not.east_exchng) then
        imax=iend
      else
        imax=iend+1
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        jmax=jend
      else
        jmax=jend+1
      endif
      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1._8/3._8)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)
      call alfabeta_tile (istr,iend,jstr,jend, imin,imax,
     &                             jmin,jmax, Bosol,Bo)
      do j=jmin,jmax
        do i=imin,imax
          Bo(i,j)=g*( Bosol(i,j)*(stflx(i,j,itemp)-srflx(i,j))
     &                              -Bo(i,j)*stflx(i,j,isalt)
     &                                                        )
          Bosol(i,j)=g*Bosol(i,j)*srflx(i,j)
          ustar(i,j)=sqrt( sqrt( 0.333333333333_8*(
     &      sustr(i,j)**2 + sustr(i+1,j)**2 + sustr(i,j)*sustr(i+1,j)
     &     +svstr(i,j)**2 + svstr(i,j+1)**2 + svstr(i,j)*svstr(i,j+1)
     &                                                            )))
          hbl(i,j)=hbls(i,j,nstp)
          kbls(i)=0
          Cr(i,N)=0._8
          Cr(i,0)=0._8
          FC(i,N)=0._8
        enddo
        do k=N-1,1,-1
          do i=imin,imax
            cff_up=(z_w(i,j,N)-z_w(i,j,k))**2
            Kern=cff_up/(cff_up +(epssfc*hbl(i,j))**2)
            FC(i,k)=FC(i,k+1) + Kern*(
     &                0.5_8*( ( u(i,j,k+1,nstp)+u(i+1,j,k+1,nstp)
     &                       -u(i,j,k  ,nstp)-u(i+1,j,k  ,nstp) )**2
     &                     +( v(i,j,k+1,nstp)+v(i,j+1,k+1,nstp)
     &                       -v(i,j,k  ,nstp)-v(i,j+1,k  ,nstp) )**2
     &                      )/(Hz(i,j,k)+Hz(i,j,k+1))
     &               -0.5_8*(Hz(i,j,k)+Hz(i,j,k+1))*( Ri_inv*bvf(i,j,k)
     &                                            +C_Ek*f(i,j)*f(i,j)
     &                                                             ))
          enddo
        enddo
        do i=imin,imax
          z_bl=z_w(i,j,0)+0.25_8*Hz(i,j,1)
          cff_up=(z_w(i,j,N)-z_bl)**2
          Kern=cff_up/(cff_up +(epssfc*hbl(i,j))**2)
            FC(i,0)=FC(i,1) + Kern*(
     &                   0.5_8*( (u(i,j,1,nstp)+u(i+1,j,1,nstp))**2
     &                        +(v(i,j,1,nstp)+v(i,j+1,1,nstp))**2
     &                       )/Hz(i,j,1)
     &                  -0.5_8*Hz(i,j,1)*( Ri_inv*bvf(i,j,1)
     &                                  +C_Ek*f(i,j)*f(i,j)
     &                                                       ))
        enddo
        do k=N,1,-1
          do i=imin,imax
            wrk1(i,k)=sqrt(swr_frac(i,j,k)*swr_frac(i,j,k-1))
            zscale=z_w(i,j,N)-z_r(i,j,k)
            Bfsfc=Bo(i,j)+Bosol(i,j)*(1._8-wrk1(i,k))
          if (Bfsfc < 0._8) zscale=min(zscale, hbl(i,j)*epssfc)
          zscale=zscale*rmask(i,j)
          zetahat=vonKar*zscale*Bfsfc
          ustar3=ustar(i,j)**3
          if (zetahat >= 0._8) then
            ws=vonKar*ustar(i,j)*ustar3/max(ustar3+5._8*zetahat, 1.D-20)
          elseif (zetahat > zeta_s*ustar3) then
            ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
          else
            ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
          endif
            Vtsq=Vtc*ws*sqrt(max(0._8, bvf(i,j,k-1) ))
            Cr(i,k)=FC(i,k)+Vtsq
            if (kbls(i) == 0 .and. Cr(i,k) < 0._8) kbls(i)=k
          enddo
        enddo
        do i=imin,imax
          if (kbls(i) > 0) then
            k=kbls(i)
            if (k == N) then
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,N)
            else
              hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*Cr(i,k+1)
     &                              -z_r(i,j,k+1)*Cr(i,k)
     &                              )/(Cr(i,k+1)-Cr(i,k))
            endif
          else
            hbl(i,j)=z_w(i,j,N)-z_w(i,j,0)
          endif
          hbl(i,j)=hbl(i,j)*rmask(i,j)
        enddo
      enddo
      cff=1.D0/12.D0 ; cff1=3.D0/16.D0
      if (istr==iwest .and. .not.west_exchng) then
        do j=jmin,jmax
          hbl(istr-1,j)=hbl(istr,j)
        enddo
      endif
      if (iend==ieast .and. .not.east_exchng) then
        do j=jmin,jmax
          hbl(iend+1,j)=hbl(iend,j)
        enddo
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        do i=imin,imax
          hbl(i,jstr-1)=hbl(i,jstr)
        enddo
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        do i=imin,imax
          hbl(i,jend+1)=hbl(i,jend)
        enddo
      endif
      if (istr==iwest .and. .not.west_exchng .and. jstr==jsouth .and. 
     &                      .not.south_exchng) then
        hbl(istr-1,jstr-1)=hbl(istr,jstr)
      endif
      if (istr==iwest .and. .not.west_exchng .and. jend==jnorth .and. 
     &                      .not.north_exchng) then
        hbl(istr-1,jend+1)=hbl(istr,jend)
      endif
      if (iend==ieast .and. .not.east_exchng .and. jstr==jsouth .and. 
     &                      .not.south_exchng) then
        hbl(iend+1,jstr-1)=hbl(iend,jstr)
      endif
      if (iend==ieast .and. .not.east_exchng .and. jend==jnorth .and. 
     &                      .not.north_exchng) then
        hbl(iend+1,jend+1)=hbl(iend,jend)
      endif
      do j=jstr-1,jend+1
        do i=istr,iend+1
          FX(i,j)=(hbl(i,j)-hbl(i-1,j))
     &                      *umask(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1
          FE1(i,j)=(hbl(i,j)-hbl(i,j-1))
     &                      *vmask(i,j)
        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)
     &                            -FX(i  ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )
     &                           -FE1(i,j  )-FE1(i-1,j+1))
        enddo
        do i=istr,iend
          hbl(i,j)=hbl(i,j) + cff1*( FX(i+1,j)-FX(i,j)
     &                              +FE(i,j+1)-FE(i,j))
          hbl(i,j)=hbl(i,j)*rmask(i,j)
        enddo
      enddo
      do j=jstr,jend
        if (.not. iic==forw_start) then
          do i=istr,iend
          enddo
        endif
        do i=istr,iend
          kbls(i)=N
        enddo
        do k=N-1,1,-1
          do i=istr,iend
            if (z_w(i,j,k) > z_w(i,j,N)-hbl(i,j)) kbls(i)=k
          enddo
        enddo
        do i=istr,iend
          k=kbls(i)
          z_bl=z_w(i,j,N)-hbl(i,j)
          zscale=hbl(i,j)
          if (swr_frac(i,j,k-1) > 0._8) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1._8 -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                     ) )
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif
            if (Bfsfc<0._8) zscale=min(zscale, hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat >= 0._8) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5._8*zetahat,
     &                                                   1.D-20 )
              ws=wm
            else
              if (zetahat > zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16._8*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat > zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
          f1=5.0_8 * max(0._8, Bfsfc) * vonKar/(ustar(i,j)**4+EPS)
          cff=1._8/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)
          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbl(i,j)*wm+EPS)
          dGm1dS(i)=min(0._8, Av_bl*f1-dAv_bl/(wm+EPS))
          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbl(i,j)*ws+EPS)
          dGt1dS(i)=min(0._8, At_bl*f1-dAt_bl/(ws+EPS))
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbl(i,j)*ws+EPS)
          dGs1dS(i)=min(0._8, As_bl*f1-dAs_bl/(ws+EPS))
          Bfsfc_bl(i)=Bfsfc
        enddo
        do i=istr,iend
          do k=N,0,-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,N)-z_w(i,j,k)
            if (Bfsfc<0._8) zscale=min(zscale, hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat >= 0._8) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5._8*zetahat,
     &                                                   1.D-20 )
              ws=wm
            else
              if (zetahat > zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16._8*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat > zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
            ssgm=(z_w(i,j,N)-z_w(i,j,k))/max(hbl(i,j),EPS)
            if (ssgm < 1._8) then
              if (ssgm<0.07D0) then
                cff=0.5_8*(ssgm-0.07D0)**2/0.07D0
              else
                cff=0.D0
              endif
              cff=cff + ssgm*(1._8-ssgm)**2
              Kv(i,j,k)=Kv(i,j,k) + wm*hbl(i,j)*cff
              Kt(i,j,k)=Kt(i,j,k) + ws*hbl(i,j)*cff
              Ks(i,j,k)=Ks(i,j,k) + ws*hbl(i,j)*cff
              if (Bfsfc < 0._8) then
                ghat(i,j,k)=Cg * ssgm*(1._8-ssgm)**2
              else
                ghat(i,j,k)=0._8
              endif
            else
              ghat(i,j,k)=0._8
              if (bvf(i,j,k) < 0._8) then
                Kv(i,j,k)=Kv(i,j,k) + nu0c
                Kt(i,j,k)=Kt(i,j,k) + nu0c
                Ks(i,j,k)=Ks(i,j,k) + nu0c
              endif
            endif
          enddo
        enddo
        do i=istr,iend
          if (rmask(i,j) > 0.5_8) then
            if (iic==forw_start) then
              do k=0,N
                Akv(i,j,k)=Kv(i,j,k)
                Akt(i,j,k,itemp)=Kt(i,j,k)
                Akt(i,j,k,isalt)=Ks(i,j,k)
              enddo
            else
              do k=0,N
                Akv(i,j,k)       = 0.5_8*Akv(i,j,k)       + 
     &                          0.5_8*Kv(i,j,k)
                Akt(i,j,k,itemp) = 0.5_8*Akt(i,j,k,itemp) + 
     &                          0.5_8*Kt(i,j,k)
                Akt(i,j,k,isalt) = 0.5_8*Akt(i,j,k,isalt) + 
     &                          0.5_8*Ks(i,j,k)
              enddo
            endif
          else
            do k=0,N
              Akv(i,j,k)=0._8
              Akt(i,j,k,itemp)=0._8
              Akt(i,j,k,isalt)=0._8
            enddo
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          hbls(i,j,3-nstp)=hbl(i,j)
        enddo
      enddo
      if (istr==iwest .and. .not.west_exchng) then
        do j=jstr,jend
          hbls(istr-1,j,3-nstp)=hbls(istr,j,3-nstp)
        enddo
      endif
      if (iend==ieast .and. .not.east_exchng) then
        do j=jstr,jend
          hbls(iend+1,j,3-nstp)=hbls(iend,j,3-nstp)
        enddo
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        do i=istr,iend
          hbls(i,jstr-1,3-nstp)=hbls(i,jstr,3-nstp)
        enddo
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        do i=istr,iend
          hbls(i,jend+1,3-nstp)=hbls(i,jend,3-nstp)
        enddo
      endif
      if (istr==iwest .and. .not.west_exchng .and. jstr==jsouth .and. 
     &                      .not.south_exchng) then
        hbls(istr-1,jstr-1,3-nstp)=hbls(istr,jstr,3-nstp)
      endif
      if (istr==iwest .and. .not.west_exchng .and. jend==jnorth .and. 
     &                      .not.north_exchng) then
        hbls(istr-1,jend+1,3-nstp)=hbls(istr,jend,3-nstp)
      endif
      if (iend==ieast .and. .not.east_exchng .and. jstr==jsouth .and. 
     &                      .not.south_exchng) then
        hbls(iend+1,jstr-1,3-nstp)=hbls(iend,jstr,3-nstp)
      endif
      if (iend==ieast .and. .not.east_exchng .and. jend==jnorth .and. 
     &                      .not.north_exchng) then
        hbls(iend+1,jend+1,3-nstp)=hbls(iend,jend,3-nstp)
      endif
      call exchange_2_tile( istr,iend,jstr,jend, Akv, N+1,
     &                      hbls(-1,-1,3-nstp),1)
      call exchange_2_tile( istr,iend,jstr,jend,
     &                   Akt(-1,-1,0,itemp), N+1,
     &                   Akt(-1,-1,0,isalt), N+1)
      end
      subroutine check_kpp_switches (ierr)
      implicit none
      integer(kind=4) ierr, is,ie, lenstr
      integer(kind=4), parameter ::
     &              LLm=1504,  MMm=1024 ,N=50
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , NT=2
      integer(kind=4), parameter :: max_opt_size=2048
      character(len=max_opt_size) cpps, srcs, kwds
      common /strings/ cpps, srcs, kwds
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+10
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<lmd_kpp.F>'
      is=ie+2 ; ie=is+16
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='INT_AT_RHO_POINTS'
      is=ie+2 ; ie=is+9
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='SMOOTH_HBL'
      is=ie+2 ; ie=is+18
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='LIMIT_UNSTABLE_ONLY'
      return
  99  if (mynode==0) write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
