clc; clear

%%
%  Generates the i,j locations of data extraction objects
%  for ROMS simulation.
%
%  Edit and rerun to append additional moorings to the same file.
%
%  Make sure the output file points to the same file you used for online
%  child boundary generation if you are also using that.
%
%  Writes to netcdf file the i and j locations of the places
%  from where we want to save data. Index locations are in [0,nx], [0,ny]
%
%  It writes a data extraction object for each set of i,j points; each boundary
%  has 3 different objects, for rho,u, and v-points. The velocity objects
%  also include and angle to which the desired velocties will be rotated
%
%  If the child grid point is not in the parent domain it is given a value
%  of -1e5
%
%  Note the mod statements for lonc and lonp. This is an attempt to deal
%  with parent and child longitudes that are possibly 360 apart
%  It will fail for grids that straddle the dateline of the zero meridian
%  In that case, subtract 180 degrees first
%

%%
% -- START USER INPUT ----------
% Parent grid directory and file name
pdir    = '../../../../Examples/USWC_model/input_data/';
pname   = 'sample_grd_riv.nc';

% Mooring details
lon_actual = -120.7;
lat_actual =  34.6;
% lat_actual =  34.56;
lon = lon_actual + 360;  % adjust withing 0 to 360 if required
lat = lat_actual + 0;    % adjust withing 0 to 360 if required

% gname = 'mooring1';
gname = 'mooring2';
% period = 40;
period = 80;
ang = 0;
mooring_vars = 'zeta, temp, salt, u, v' ;

% Output file name and info
ename   = 'sample_edata.nc';
info    = ['indices for ' gname ' in ' pname];

% -- END USER INPUT ------------
%%

pname = [pdir pname];

lonp = ncread(pname,'lon_rho');
latp = ncread(pname,'lat_rho');
lonp = mod(lonp,360);

obj_name = gname;
obj_lon = lon;
obj_lat = lat;
obj_ang = ang;
obj_msk =   1;
add_object(ename,obj_name,gname,lonp,latp,period,obj_lon,obj_lat,obj_msk,obj_ang);
ncwriteatt(ename,obj_name,'output_vars',mooring_vars);
ncwriteatt(ename,obj_name,'lat',lat_actual);
ncwriteatt(ename,obj_name,'lon',lon_actual);
ncwriteatt(ename, '/', [gname '_info'],  info);           % info on parent and child grid

