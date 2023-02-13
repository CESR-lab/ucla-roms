
	function fld = frc_era(data,grd,varname,irec,method)

	datname  = data.datname;
	lon_frc = data.lon;
	lat_frc = data.lat;
	i0 = data.i0;
	j0 = data.j0;
	fnx = data.fnx;
	fny = data.fny;

	lon = grd.lon;
	lat = grd.lat;
	mask= grd.mask;

        frc = ncread(datname,varname,[i0 j0 irec],[fnx fny 1]); 
	frc = fliplr(frc);

        frc(mask<1) = nan;
        frc = inpaint_nans(frc,2);
        fld = interp2(lon_frc,lat_frc,frc',lon,lat,method);

	return
