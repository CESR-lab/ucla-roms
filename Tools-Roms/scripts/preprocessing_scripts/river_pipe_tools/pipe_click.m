
%load /avatar/nmolem/batavia/nmolem/OBSERV/COAST/coast_pac

 lon  = ncread(grdname,'lon_rho')';
 lat  = ncread(grdname,'lat_rho')';
 mask = ncread(grdname,'mask_rho')';
 rflx = ncread(grdname,'pipe_flux')';

 lon = lon(j0:j1,i0:i1);
 lat = lat(j0:j1,i0:i1);
 mask= mask(j0:j1,i0:i1);
 rflx= rflx(j0:j1,i0:i1);

 lonp = 0.25*(lon(1:end-1,1:end-1)+lon(2:end,1:end-1)+lon(1:end-1,2:end)+lon(2:end,2:end));
 latp = 0.25*(lat(1:end-1,1:end-1)+lat(2:end,1:end-1)+lat(1:end-1,2:end)+lat(2:end,2:end));

 sc1 = ipip + 1;
 sc0 =-2*sc1/255;

 i = 1;j=1;
 h = figure('WindowButtonDownFcn',{@river_wbdcb,i,j,lon,lat,lonp,latp,rflx,mask,sc0,sc1});

nmask = 0*mask;
nmask(~mask) = -1e5;
nmask(rflx>0) = 0;
mypcolor(lon(2:end-1,2:end-1),lat(2:end-1,2:end-1),rflx(2:end-1,2:end-1)+nmask(2:end-1,2:end-1))
clear jet
cm = colormap(jet(256));
cm(1,:) = [204 153 0]/255;
colormap(cm);
caxis([sc0 sc1]);
hold on
plot(pip_lon,pip_lat,'g*')
title('Close figure when done with this pipe')

 uiwait(h);
 disp(['Done with pipe: ' num2str(ipip)])

%close all
%imagesc(rflx);axis xy
 ncwrite(grdname,'pipe_flux',rflx',[i0 j0]);

