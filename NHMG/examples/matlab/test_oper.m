%% Test elliptic operator

sizex = 5e3*2;
sizez = 2e3*2;

resmx = 3;
nzmx = 40*2^(resmx-1);
prof_z = zeros(nzmx,resmx);
prof_u = zeros(nzmx,resmx); 
Ekconv = zeros(resmx,1);

for res = resmx:resmx %% test resolution convergence
nx = 100*2^(res-1);
nz =  40*2^(res-1);

dx = sizex/nx; dxi = 1./dx; dxi2 = dxi^2;
dz = sizez/nz; dzi = 1./dz; dzi2 = dzi^2;

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
mask = ones(nx,nz);
masku = zeros(nx+1,nz);
maskw = zeros(nx,nz+1);
rhs = zeros(nx,nz);
for k = 1:nz
   for i=1:nx
      mask(i,k) = z1(k)> h1(i);
   end
end
maskw(:,2:nz) = mask(:,1:nz-1)&mask(:,2:nz);maskw(:,nz+1) = mask(:,nz);
masku(2:nx,:) = mask(1:nx-1,:)&mask(2:nx,:);

bet = 600/sizex^2;
x0 = sizex*0.65;
z0 = sizez*(0.75-1);
r1 = (xr-x0).^2 + (zr-z0).^2;
x0 = sizex*0.75;
z0 = sizez*(0.65-1);
r2 = (xr-x0).^2 + (zr-z0).^2;
rhs = exp(-bet*r1) - exp(-bet*r2);
rhs(mask<1) = 0.0;

u = zeros(nx+1,nz);
w = cos(xw*pi/sizex); w = w.*maskw;
rhs = (u(2:nx+1,:)-u(1:nx,:))*dxi + (w(:,2:nz+1)-w(:,1:nz))*dzi;
div1 =rhs;

w0 = w;
u0 = u;

Ek0 = sum(sum(w(:,2:nz).^2)) + sum(0.5*w(:,nz+1).^2) + sum(sum(u.^2));
Ek0 = 0.5*Ek0*dz*dx/(sizex*sizez)

ny = 1; dy = dx; sizey = dy;
u0 = zeros(nx+1,ny,nz); u0(:,1,:) = u;
v0 = zeros(nx,ny+1,nz); v = v0;
w0 = zeros(nx,ny,nz+1); w0(:,1,:) = w;
ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);
nd = ndu+ndv+ndw;

X = zeros(nd,1);
X(1:ndu) = u0(:);
X(ndu+1:ndu+ndv) = v0(:);
X(ndu+ndv+1:end) = w0(:);

Au = ones(nx+1,ny,nz)/dx;
Av = ones(nx,ny+1,nz)/dy;
Aw = ones(nx,ny)/dz;

Mzl = makeMzl(dx,dy,dz,nx,ny,nz);
mask3 = zeros(nx,ny,nz);
mask3(:,1,:) = mask;
Dzl = makeDzl(dx,dy,dz,mask3);
A2 = Dzl*Dzl';

Ek0 = 0.5*X'*Mzl*X/(sizex*sizey*sizez)

nd = nx*nz;
tic
if 0
A = sparse(nd,nd,5*nd);
is = 1;
ks = nx;
for i = 1:nx
   for k = 1:nz
      idx = 1*((k-1)*nx + i-1);
      diag = 0;
      if mask(i,k)
       if i>1
        if mask(i-1,k)
          diag = diag - dxi2;
          A(idx+1,idx+1-is) = dxi2;
        end
       end
       if i<nx
        if mask(i+1,k)
          diag = diag - dxi2;
          A(idx+1,idx+1+is) = dxi2;
        end
       end
       if k>1
        if mask(i,k-1)
          diag = diag - dzi2;
          A(idx+1,idx+1-ks) = dzi2;
        end
       end
       if k<nz
        if mask(i,k+1)
          diag = diag - dzi2;
          A(idx+1,idx+1+ks) = dzi2;
        end
       else
          diag = diag-2*dzi2;
       end
      end
      if diag==0
         A(idx+1,idx+1) = 1;
      else
         A(idx+1,idx+1) = diag;
      end
   end
end
else

ce = zeros(nx,nz); ce(1:nx-1,:) = dxi2;
cw = zeros(nx,nz); cw(2:nx  ,:) = dxi2;
ct = zeros(nx,nz); ct(:,1:nz-1) = dzi2;
cb = zeros(nx,nz); cb(:,2:nz  ) = dzi2;
cw(2:nx  ,:) = mask(1:nx-1,:).*cw(2:nx  ,:).*mask(2:nx  ,:);
ce(1:nx-1,:) = mask(2:nx  ,:).*ce(1:nx-1,:).*mask(1:nx-1,:);
cb(:,2:nz  ) = mask(:,1:nz-1).*cb(:,2:nz  ).*mask(:,2:nz  );
ct(:,1:nz-1) = mask(:,2:nz  ).*ct(:,1:nz-1).*mask(:,1:nz-1);
cc = -ce-cw-ct-cb;
cc(:,nz) = cc(:,nz) - 2*dzi2;
cc(mask<1) = -1.0;

id = [-nx -1 0 1 nx];
diags = zeros(nd,5);
cb = reshape(cb,[nd,1]);diags(1:end-nx,1) = cb(nx+1:end);
cw = reshape(cw,[nd,1]);diags(1:end-1 ,2) = cw(1+1:end);
ce = reshape(ce,[nd,1]);diags(2:end   ,4) = ce(1:end-1);
ct = reshape(ct,[nd,1]);diags(nx+1:end,5) = ct(1:end-nx);
diags(:,3) = reshape(cc,[nd,1]);
A = spdiags(diags,id,nd,nd);
end
toc

tic
rhs2 = reshape(rhs,[nx*nz,1]);
p_b = A\rhs2;
p2 = reshape(p_b,nx,nz);
toc



u(2:nx,:) = u(2:nx,:) - (p2(2:nx,:) - p2(1:nx-1,:))*dxi;
w(:,2:nz) = w(:,2:nz) - (p2(:,2:nz) - p2(:,1:nz-1))*dzi;
w(:,nz+1) = w(:,nz+1) + 2*p2(:,nz)*dzi;
w = w.*maskw;
u = u.*masku;
div = (u(2:nx+1,:)-u(1:nx,:))*dxi + (w(:,2:nz+1)-w(:,1:nz))*dzi;
Ek2 = sum(sum(w(:,2:nz).^2)) + sum(0.5*w(:,nz+1).^2) + sum(sum(u.^2));
Ek2 = 0.5*Ek2*dz*dx/(sizex*sizez)

X(1:ndu) = u(:);
X(ndu+1:ndu+ndv) = v(:);
X(ndu+ndv+1:end) = w(:);
Ek2 = 0.5*X'*Mzl*X/(sizex*sizey*sizez)

   prof_z(1:nz,res) = zu(nx/2+1,:);
   prof_u(1:nz,res) =  u(nx/2+1,:); 
   Ekconv(res) = Ek2;
end %% convergence loop
[max(p_b) min(p_b)];
p2(mask<1) = nan;

figure(2)
pcolor(xw,zw,w);shading flat
title('z coord')
colorbar
hold on
plot(x1,h1,'k','linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',15)
wmx = max(max(w));
contour(xw,zw,w,[-1:0.1:1]*wmx,'k')
hold off
