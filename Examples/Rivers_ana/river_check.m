%% Checking integrals of water input and T, and S conservation

fname = 'ptest_his.0000.nc';
fname = 'uswc_his.0000.nc';

nx = 100;
ny = 100;
nz = 10;

h = ncread(fname,'h',[2,2],[nx,ny]);
mask = ncread(fname,'mask_rho',[2,2],[nx,ny]);

ptim = ncread(fname,'ocean_time');
nt = length(ptim);

dxi = ncread(fname,'pm');
dyi = ncread(fname,'pn');

dx = 1./mean(dxi(:));
dy = 1./mean(dyi(:));

ts = ncreadatt(fname,'/','theta_s');
tb = ncreadatt(fname,'/','theta_b');
hc = ncreadatt(fname,'/','hc');

vol = zeros(nt,1);
for it = 1:nt
  zeta = ncread(fname,'zeta',[2,2,it],[nx,ny,1]);
  zeta(mask<1) = nan;
  vol(it) = dx*dy*nansum(zeta(:));
end

s_int = zeros(nt,1);
t_int = zeros(nt,1);

for it = 1:nt
  zeta = ncread(fname,'zeta',[2,2,it],[nx,ny,1]);
  salt = ncread(fname,'salt',[2,2,1,it],[nx,ny,nz,1]);
  temp = ncread(fname,'temp',[2,2,1,it],[nx,ny,nz,1]);
% temp = temp-10;
% salt = salt-36;

  zw = zlevs3(h,zeta,ts,tb,hc,nz,'w','new2008');
  zw = permute(zw,[2 3 1]);
  dz = zw(:,:,2:end)-zw(:,:,1:end-1);
  t_int(it) = dx*dy*nansum(temp(:).*dz(:));
  s_int(it) = dx*dy*nansum(salt(:).*dz(:));
end

dt = gradient(t_int,600);
ds = gradient(s_int,600);
dv = gradient(vol,600);
