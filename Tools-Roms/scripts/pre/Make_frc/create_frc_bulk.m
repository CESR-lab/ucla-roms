function create_frc_bulk(gridfile,frcfile)
%
%   Create ROMS bulk forcing file
% 
%
[nx,ny] = size(ncread(gridfile,'h'));

%
%  Create variables
%
nccreate(frcfile,'time','Dimensions',{'time',0},'datatype','double');
ncwriteatt(frcfile,'time','long_name','time since 2000-1-1');
ncwriteatt(frcfile,'time','units','day');

nccreate(frcfile,'uwnd','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'uwnd','long_name','10 meter wind in x-direction');
ncwriteatt(frcfile,'uwnd','units','m/s');

nccreate(frcfile,'vwnd','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'vwnd','long_name','10 meter wind in y-direction');
ncwriteatt(frcfile,'vwnd','units','m/s');

nccreate(frcfile,'Tair','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'Tair','long_name','air temperature at 2m');
ncwriteatt(frcfile,'Tair','units','degrees C');

nccreate(frcfile,'swrad','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'swrad','long_name','downward solar shortwave radiation');
ncwriteatt(frcfile,'swrad','units','W/m^2');

nccreate(frcfile,'lwrad','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'lwrad','long_name','downward long wave radiation');
ncwriteatt(frcfile,'lwrad','units','W/m^2');

nccreate(frcfile,'qair','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'qair','long_name','Relative humidity at 2 m');
ncwriteatt(frcfile,'qair','units','kg/kg');

nccreate(frcfile,'rain','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',0},'datatype','single');
ncwriteatt(frcfile,'rain','long_name','Total precipitation');
ncwriteatt(frcfile,'rain','units','m/s');


%
%
%  Write global attributes
%
 ncwriteatt(frcfile,'/','Title','ROMS Bulk surface forcing file');
 ncwriteatt(frcfile,'/','Date',date);
 ncwriteatt(frcfile,'/','gridfile',gridfile);
%
%
return


