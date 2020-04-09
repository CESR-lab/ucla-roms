%% Test zlevel operator

sizex = 1e4;
sizey = 1e2;
sizez = 4e3; 

resmx = 7;
nzmx = 40*2^(resmx-1);
prof_z = zeros(nzmx,resmx);
prof_u = zeros(nzmx,resmx); 
Ekconv = zeros(resmx,1);

for res = 5:resmx %% test resolution convergence
  nx = 50*2^(res-1);
  ny = 1;
  nz = 20*2^(res-1);
  
  ndu = (nx+1)*ny*nz;
  ndv = nx*(ny+1)*nz;
  ndw = nx*ny*(nz+1);
  nd = ndu+ndv+ndw;
  
  dx = sizex/nx; 
  dy = sizey/ny;
  dz = sizez/nz;

  x1 = sizex*[0.5:1:nx-0.5]/nx;
  z1 = sizez*[0.5:1:nz-0.5]/nz - sizez;
  x1u = sizex*[0.0:1:nx-0.0]/nx;
  z1w = sizez*[0:nz]/nz - sizez;
  [xr,zr] = meshgrid(x1,z1);
  [xu,zu] = meshgrid(x1u,z1);
  [xw,zw] = meshgrid(x1,z1w);
  xr = xr';
  zr = zr';
  xu = xu';
  zu = zu';
  xw = xw';
  zw = zw';

  alpha = 25./(sizex^2);
  h1 = -(sizez - 0.5*sizez*exp(-alpha*(x1-0.5*sizex).^2));
  mask = ones(nx,ny,nz);

  for k = 1:nz
    mask(:,1,k) = h1<z1(k);
  end
  
  Mzl = makeMzl(dx,dy,dz,nx,ny,nz);
  ndt = ndu+ndv+ndw;
  Mizl = spdiags(1./spdiags(Mzl),0,ndt,ndt);
  Dzl = makeDzl(dx,dy,dz,mask);
  
  Pzl = Dzl*Mizl*Dzl';

  u = zeros(nx+1,ny,nz);
  v = zeros(nx,ny+1,nz);
  w = zeros(nx,ny,nz+1);
  w(:,1,:) = cos(xw*pi/sizex); 
  
  w(:,:,1) = 0;
  w(:,:,2:nz) = w(:,:,2:nz).*mask(:,:,1:nz-1).*mask(:,:,2:nz);
  w(:,:,nz+1) = w(:,:,nz+1).*mask(:,:,nz);
  
  X = zeros(nd,1);
  X(1:ndu) = u(:);
  X(ndu+1:ndu+ndv) = v(:);
  X(ndu+ndv+1:end) = w(:);
  
  R = Dzl*X;
  p = Pzl\R;
  dp = Mizl*Dzl'*p(:);

  Xd = X-dp; 

  Ek = 0.5*Xd'*Mzl*Xd/(sizex*sizey*sizez)

  ud = reshape(Xd(1:ndu),nx+1,ny,nz);
  wd = reshape(Xd(ndu+ndv+1:end),nx,ny,nz+1);
  Ekconv(res) = Ek;
  
  load ref_grid
  wi = interp2(xw',zw',squeeze(wd)',xw_ref,zw_ref);
  wi = squeeze(wi);
  fname = ['Zlev' num2str(res-1)];
  save(fname,'wi');
end

wd = squeeze(wd);
figure(2)
pcolor(xw,zw,wd);shading flat
title('z coord')
colorbar
hold on
plot(x1,h1,'k','linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',15)
wmx = max(max(wd));
contour(xw,zw,wd,[-1:0.1:1]*wmx,'k')
hold off
