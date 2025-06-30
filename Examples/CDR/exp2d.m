dx = 1e-3;

lon = [-10:dx:10];
lat = lon;
[lon2,lat2] = meshgrid(lon,lat);
lon2 = lon2';
lat2 = lat2';

d2 = (lon2.^2+ lat.^2);
scl2 = 3*3;

func = exp(-d2/scl2);

imagesc(func')
colorbar

sum(func(:))*dx*dx
