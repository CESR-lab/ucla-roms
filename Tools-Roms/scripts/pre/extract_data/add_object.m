function add_object(ename,obj_name,gname,lonp,latp,period,obj_lon,obj_lat,obj_mask,obj_ang)
  % Adds a data extraction object to a netcdf file
  % The object contains the fractional i,j positions of the desired locations
  % plus optionally an angle with the desired rotation for velocities

  disp(['Adding: ', obj_name, ' to ',ename]);

  [nxp,nyp] = size(lonp);                  % dimension sizes of parent domain
% ip_rho = [0:nxp-1]-0.5;                  % index numbering of parent domain
% jp_rho = [0:nyp-1]-0.5;
% [ip_rho,jp_rho] = meshgrid(ip_rho,jp_rho);
% ip_rho = ip_rho';
% jp_rho = jp_rho';


  % crop parent grid to minimal size
  [lons,lats,ips,jps] = crop_parent(lonp,latp,obj_lon,obj_lat);

  % interpolate the parent indices:
  obj_i = griddata(lons,lats,ips,obj_lon,obj_lat);
  obj_j = griddata(lons,lats,jps,obj_lon,obj_lat);

  obj_i(obj_mask<1) = -1e5;
  obj_j(obj_mask<1) = -1e5;

  if sum(isnan(obj_i))>0
    disp(obj_name)
    disp('Some points are outside the grid and not masked!!')
    error 'fatal'
  end
  i_chk = obj_i; i_chk(obj_mask<1) = [];
  if (min(i_chk)<0)|(max(i_chk)>nxp-2)
    disp('Some points are borderline')
  end
  j_chk = obj_j; j_chk(obj_mask<1) = [];
  if (min(j_chk)<0)|(max(j_chk)>nyp-2)
    disp('Some points are borderline')
  end
              
  obj_i(isnan(obj_i)) = -1e5;  
  obj_j(isnan(obj_j)) = -1e5;  

  np = length(obj_i);

  % code to help with nested boundaries objects
  if strfind(obj_name,'grid')
    if ~isempty(strfind(obj_name,'west')) |...
       ~isempty(strfind(obj_name,'east'))
      if strcmp(obj_name(end),'v')
        dname = [gname '_eta_v'];
      else
        dname = [gname '_eta_rho'];
      end
    elseif ~isempty(strfind(obj_name,'north')) |...
           ~isempty(strfind(obj_name,'south'))
      if strcmp(obj_name(end),'u')
        dname = [gname '_xi_u'];
      else
        dname = [gname '_xi_rho'];
      end
    else
      disp(obj_name)
      error 'Grid bry object, but not an obvious boundary'
    end
  else
   dname = [gname '_np'];
  end

  if exist('obj_ang')
    nccreate(ename,obj_name, 'dimensions', {dname,np,'three',3}, 'datatype', 'single');
    ncwriteatt(ename,obj_name, 'long_name', 'index coordinates and angle of data object');
    ncwriteatt(ename,obj_name, 'units', 'non-dimensional and radian');
    ncwrite(ename,obj_name, [obj_i obj_j obj_ang]);
  else
    nccreate(ename,obj_name, 'dimensions', {dname,np,'two',2}, 'datatype', 'single');
    ncwriteatt(ename,obj_name, 'long_name', 'index coordinates of data object');
    ncwriteatt(ename,obj_name, 'units', 'non-dimensional');
    ncwrite(ename,obj_name, [obj_i obj_j]);
  end

  ncwriteatt(ename,obj_name,'output_period',period);


  return


