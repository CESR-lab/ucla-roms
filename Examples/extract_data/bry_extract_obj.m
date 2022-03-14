%  Generates the i,j locations of data extraction objects
%  for use as boundary forcing in a subsequent nested grid
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

% -- START USER INPUT ----------
% Parent grid directory and file name
pdir    = './';
pname   = 'sample_grd_river.nc';

% Child grid directory and file name
cdir    = './';
cname   = 'child_sample_grd.nc';

% Output file name and info
ename   = 'sample_edata.nc';
info    = ['indices for ' cname ' in ' pname];

% Open boundaries: south,west,north,east. 1=yes, 0=no.
OBC_SWNE = [1,1,1,0];

 gname = 'grid1'
 obc_west = 1;
 obc_east = 0;
 obc_south= 1;
 obc_north= 1;
 period = 1800;


% -- END USER INPUT ------------

delete(ename);                            % prevents overwriting error

pname = [pdir pname];
cname = [cdir cname];

lonp = ncread(pname,'lon_rho');        
latp = ncread(pname,'lat_rho');
lonp = mod(lonp,360);

lonr = ncread(cname,'lon_rho');
latr = ncread(cname,'lat_rho');
angr = ncread(cname,'angle');
lonr = mod(lonr,360);

lonu = 0.5*(lonr(1:end-1,:) + lonr(2:end,:));
latu = 0.5*(latr(1:end-1,:) + latr(2:end,:));
angu = 0.5*(angr(1:end-1,:) + angr(2:end,:));
lonv = 0.5*(lonr(:,1:end-1) + lonr(:,2:end));
latv = 0.5*(latr(:,1:end-1) + latr(:,2:end));
angv = 0.5*(angr(:,1:end-1) + angr(:,2:end));


  if obc_west
    bnd = 'west';
    obj_name = [gname '_' bnd '_r'];
    obj_lon = lonr(1,:)';
    obj_lat = latr(1,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat);

    obj_name = [gname '_' bnd '_u'];
    obj_lon = lonu(1,:)';
    obj_lat = latu(1,:)';
    obj_ang = angu(1,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang);

    obj_name = [gname '_' bnd '_v'];
    obj_lon = lonv(1,:)';
    obj_lat = latv(1,:)';
    obj_ang = angv(1,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang);
  end

  if obc_east
    bnd = 'east';
    obj_name = [gname '_' bnd '_r'];
    obj_lon = lonr(end,:)';
    obj_lat = latr(end,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat)

    obj_name = [gname '_' bnd '_u'];
    obj_lon = lonu(end,:)';
    obj_lat = latu(end,:)';
    obj_ang = angu(end,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)

    obj_name = [gname '_' bnd '_v'];
    obj_lon = lonv(end,:)';
    obj_lat = latv(end,:)';
    obj_ang = angv(end,:)';
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)
  end

  if obc_south
    bnd = 'south';
    obj_name = [gname '_' bnd '_r'];
    obj_lon = lonr(:,1);
    obj_lat = latr(:,1);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat)

    obj_name = [gname '_' bnd '_u'];
    obj_lon = lonu(:,1);
    obj_lat = latu(:,1);
    obj_ang = angu(:,1);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)

    obj_name = [gname '_' bnd '_v'];
    obj_lon = lonv(:,1);
    obj_lat = latv(:,1);
    obj_ang = angv(:,1);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)
  end

  if obc_north
    bnd = 'north';
    obj_name = [gname '_' bnd '_r'];
    obj_lon = lonr(:,end);
    obj_lat = latr(:,end);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat)

    obj_name = [gname '_' bnd '_u'];
    obj_lon = lonu(:,end);
    obj_lat = latu(:,end);
    obj_ang = angu(:,end);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)

    obj_name = [gname '_' bnd '_v'];
    obj_lon = lonv(:,end);
    obj_lat = latv(:,end);
    obj_ang = angv(:,end);
    add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_ang)
  end

  ncwriteatt(ename, '/', 'info',  info);           % info on parent and child grid

