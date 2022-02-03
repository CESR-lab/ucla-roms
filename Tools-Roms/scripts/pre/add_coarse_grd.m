function add_coarse_grd(gname);

lon = ncread(gname,'lon_rho');
lat = ncread(gname,'lat_rho');
h   = ncread(gname,'h');
ang = ncread(gname,'angle');
mask= ncread(gname,'mask_rho');

 lonc = f2c(lon);
 latc = f2c(lat);
 hc   = f2c(h);
 angc = f2c(ang);
 maskc= f2c(mask);
 maskc(maskc>0.5)=1;
 maskc(maskc<0.5)=0;

 [nxc,nyc] = size(lonc);

 % Check if the coarse variables already exist in the grid file
 % Otherwise, create the coarse variables
 info = ncinfo(gname);
 nvars = length(info.Variables)
 present = 0;
 for i=1:nvars
   if strcmp(info.Variables(i).Name,'lon_coarse')
     present = 1;
   end
 end

 if ~present % can't add if already added
   nccreate(gname,'lon_coarse','dimensions',{'xi_coarse',nxc,'eta_coarse',nyc},'datatype','double');
   nccreate(gname,'lat_coarse','dimensions',{'xi_coarse',nxc,'eta_coarse',nyc},'datatype','double');
   nccreate(gname,'h_coarse','dimensions',{'xi_coarse',nxc,'eta_coarse',nyc},'datatype','double');
   nccreate(gname,'angle_coarse','dimensions',{'xi_coarse',nxc,'eta_coarse',nyc},'datatype','double');
   nccreate(gname,'mask_coarse','dimensions',{'xi_coarse',nxc,'eta_coarse',nyc},'datatype','double');
 end

 ncwrite(gname,'lon_coarse',lonc);
 ncwrite(gname,'lat_coarse',latc);
 ncwrite(gname,'h_coarse',hc);
 ncwrite(gname,'angle_coarse',angc);
 ncwrite(gname,'mask_coarse',maskc);

return
%--------------------------------------------------
function  fc = f2c(f)

 fc = f2c_xdir(f);
 fc = fc';
 fc = f2c_xdir(fc);
 fc = fc';

 return

%--------------------------------------------------
function  fc = f2c_xdir(f)

 [nxp,nyp] = size(f);
 nxcp = (nxp-2)/2+2;

 fc = zeros(nxcp,nyp);
 fc(2:end-1,:) = 0.5*(f(2:2:end-1,:)+f(3:2:end-1,:));
 fc(1,:) = f(1,:) + 0.5*(f(1,:)-f(2,:));
 fc(nxcp,:) = f(nxp,:) + 0.5*(f(nxp,:)-f(nxp-1,:));

 return


