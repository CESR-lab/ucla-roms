
	function [lons,lats,ips,jps] = crop_parent(lons,lats,lonbc,latbc);

	% Find limits of child grid and 
	% crop parent grid to minimal size

        % The domain index space of the parent is [0:nx]x[0:ny]
        % rho-points are from -0.5 to nx+0.5;
        % where nx = nxp-2;

        [nxp,nyp] = size(lons);             % dimension sizes of parent domain (nxp=nx+2)
        ip_rho = [0:nxp-1]-0.5;             % index numbering of parent domain
	jp_rho = [0:nyp-1]-0.5;             % starts at -0.5, ends at nx+0.5
        [ips,jps] = meshgrid(ip_rho,jp_rho);
        ips = ips';
        jps = jps';

	lonbc_mn = min(lonbc(:));
	lonbc_mx = max(lonbc(:));
	latbc_mn = min(latbc(:));
	latbc_mx = max(latbc(:));

	for it = 1:5

          [nxs,nys] = size(lons);
	  lon_mn = min(lons');
          lon_mx = max(lons');
    	  lat_mn = min(lats);
          lat_mx = max(lats);

	  i0 = find(lon_mx<lonbc_mn,1,'last');  if isempty(i0); i0= 1 ; end
	  i1 = find(lon_mn>lonbc_mx,1,'first'); if isempty(i1); i1=nxs; end
	  j0 = find(lat_mx<latbc_mn,1,'last');  if isempty(j0); j0= 1 ; end
	  j1 = find(lat_mn>latbc_mx,1,'first'); if isempty(j1); j1=nys; end

  	  lons = lons(i0:i1,j0:j1);
	  lats = lats(i0:i1,j0:j1);
  	  ips  = ips(i0:i1,j0:j1);
	  jps  = jps(i0:i1,j0:j1);
	end

	if 0
          close all
          plot(lons,lats,'.k')
          hold on
          plot(lonbc,latbc,'.r')
          hold off
	  pause(1.0)
        end

        end
