
%clear all
rname = 'FMbc'

Ekconv = zeros(6,1);

for pw = 0:0
  fac = 2^pw;
  testM3D
  return
  uc = zeros(nx+1,ny,nz);
  vc = zeros(nx,ny+1,nz);
  wc = cos(xw*pi/sizex);%wc(:,:,1);

%% A trick to get momenta from cartesian;
%% u = uc + zx*wc; and U = u - zx*w;
Xc = [uc(:)' vc(:)' -wc(:)']'; 
Xtmp = T*Xc;
u = reshape(Xtmp(1:ndu),nx+1,ny,nz);
u(2:nx,:,1) = uc(2:nx,:,1) + 0.5*(zx(2:nx,:,1).*wc(2:nx,:,2)+zx(1:nx-1,:,1).*wc(1:nx-1,:,2));
v = reshape(Xtmp(ndu+1:ndu+ndv),nx,ny+1,nz);
w = wc;

disp('solving')
tic
if strcmp(rname(1:2),'FM')
  disp('Flux Momentum')
  X = [u(:)' v(:)' w(:)']'; 
  R = D*T*X;
  p = Pt\R;
  dp = Mi*D'*p(:);  
else
  disp('Carthesian')
  X = [uc(:)' vc(:)' wc(:)']';
  R = Dc*X;
  p = Pc\R;
  dp = Mi*Dc'*p(:);
end
toc

Xd = X-dp; 

if strcmp(rname(1:2),'FM')
  Ekd = 0.5*Xd'*M*T*Xd/(sizex*sizey*sizez)
else
  Ekd = 0.5*Xd'*M*Xd/(sizex*sizey*sizez)  
end
    Pwrk = 0.5*Xd'*(M*T)'*dp/(sizex*sizey*sizez)
%    Pwrk = 0.5*dp'*M*T*Xd/(sizex*sizey*sizez)

%   Ek0 = 0.5*X'*M*X/(sizex*sizey*sizez)

Ekconv(pw+1) = Ekd;
pz = reshape(dp(ndu+ndv+1:end),nx,nz+1);
%[pz,px] = gradient(reshape(p,nx,nz));

ud = reshape(Xd(1:ndu),nx+1,ny,nz);
wd = reshape(Xd(ndu+ndv+1:end),nx,ny,nz+1);
wd(:,:,1) = 0.5*(ud(2:nx+1,:,1)+ud(1:nx,:,1)).*zx(:,:,1)./(1+zx(:,:,1).^2);
wd = squeeze(wd);
zw = squeeze(zw);
xw = squeeze(xw);



%% interpolate w to grid at 50x20 resolution

nxref = 50;
nzref = 20;
iref = [0.5:1:nxref-0.5]/nxref;
kref = [0.0:1:nzref-0.0]/nzref;
[iref,kref]=meshgrid(iref,kref);
iref = iref';
kref = kref';

id = [0.5:1:nx-0.5]/nx;
kd = [0.0:1:nz-0.0]/nz;
[id,kd]=meshgrid(id,kd);
%id = id';
%kd = kd';
wi = interp2(id,kd,wd',iref,kref);

fname = ['Data/' rname num2str(pw)];
save(fname,'wi');
fname = ['Data/Ek_' rname];
save(fname,'Ekconv');
end


figure(1)
hold off
pcolor(xw,zw,wd);shading flat
colormap(jet(256))
colorbar
hold on
plot(xw(:,1),zw(:,1),'k','linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',15)
wmx = max(max(wd));
contour(xw,zw,wd,[-1:0.1:1]*wmx,'k')
hold off

return





q = Pt\p;
qc = Pc\p;

q2 = reshape(q,nx,nz);
qc2 = reshape(qc,nx,nz);

x2 = 0*q2;
z2 = 0*q2;

for k=1:nz
  x2(:,k) = x(:,1);
  z2(:,k) = zr(:,1,k);
end