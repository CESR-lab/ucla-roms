% Check 3-dimensional fields of kinetic energy.

gname = 'input_data/sample_grd_riv.nc';
dname = 'diags_dia.20121209143435.nc';
hname = 'diags_his.20121209133435.nc';
ename = 'diags_ekb.20121209143435.nc';

ifr = 1;

if findstr(ncreadatt(dname,'/','normalization'),'include');
  diag_norm = logical(0);
  disp('This version only works for normalized diagnostics')
  return
else
  diag_norm = logical(1);
end

[nxu,nyp,nz] = size(ncread(hname,'u',[1 1 1 ifr],[inf inf inf 1]));
nx = nxu-1;
ny = nyp-2;
nyv= nyp-1;

h = ncread(gname,'h');
dx = 1./ncread(gname,'pm');
dy = 1./ncread(gname,'pn');


t1 = ncread(hname,'ocean_time',[ifr],[1]);
u1 = ncread(hname,'u',[1 2 1 ifr],[nxu ny inf 1]);
v1 = ncread(hname,'v',[2 1 1 ifr],[nx nyv inf 1]);
u1(1,:,:) = 0; u1(end,:,:) = 0;
v1(:,1,:) = 0; v1(:,end,:) = 0;

u21 = 0.5*u1.*u1;
v21 = 0.5*v1.*v1;
e1 = 0.5*(u21(2:end,:,:)+u21(1:end-1,:,:)+v21(:,2:end,:)+v21(:,1:end-1,:) );

z1 = ncread(hname,'zeta',[1 1 ifr],[inf inf 1]);
zw1 = zlevs(h,z1,hname,'w');
dz1 = zw1(:,:,2:end)-zw1(:,:,1:end-1);
dzu1 = 0.5*(dz1(2:end,2:end-1,:) + dz1(1:end-1,2:end-1,:));
dzv1 = 0.5*(dz1(2:end-1,2:end,:) + dz1(2:end-1,1:end-1,:));

t2 = ncread(hname,'ocean_time',[ifr+1],[1]);
u2 = ncread(hname,'u',[1 2 1 ifr+1],[nxu ny inf 1]);
v2 = ncread(hname,'v',[2 1 1 ifr+1],[nx nyv inf 1]);
u2(1,:,:) = 0; u2(end,:,:) = 0;
v2(:,1,:) = 0; v2(:,end,:) = 0;
u22 = 0.5*u2.*u2;
v22 = 0.5*v2.*v2;
e2 = 0.5*(u22(2:end,:,:)+u22(1:end-1,:,:)+v22(:,2:end,:)+v22(:,1:end-1,:) );

z2 = ncread(hname,'zeta',[1 1 ifr+1],[inf inf 1]);
zw2 = zlevs(h,z2,hname,'w');
dz2 = zw2(:,:,2:end)-zw2(:,:,1:end-1);
dzu2 = 0.5*(dz2(2:end,2:end-1,:) + dz2(1:end-1,2:end-1,:));
dzv2 = 0.5*(dz2(2:end-1,2:end,:) + dz2(2:end-1,1:end-1,:));

uhlf = 0.5*(u1+u2);
vhlf = 0.5*(v1+v2);
dzuhlf = 0.5*(dzu1+dzu2);
dzvhlf = 0.5*(dzv1+dzv2);
dzhlf = 0.5*(dz1+dz2);
dzhlf = dzhlf(2:end-1,2:end-1,:);


edia1 = ncread(ename,'ek_pgr',[2 2 1 ifr],[nx ny inf 1]);
edia2 = ncread(ename,'ek_cor',[2 2 1 ifr],[nx ny inf 1]);
edia3 = ncread(ename,'ek_adv',[2 2 1 ifr],[nx ny inf 1]);
edia4 = ncread(ename,'ek_dis',[2 2 1 ifr],[nx ny inf 1]);
edia5 = ncread(ename,'ek_hmx',[2 2 1 ifr],[nx ny inf 1]);
edia6 = ncread(ename,'ek_vmx',[2 2 1 ifr],[nx ny inf 1]);
edia7 = ncread(ename,'ek_cpl',[2 2 1 ifr],[nx ny inf 1]);
edia8 = ncread(ename,'ek_dzt',[2 2 1 ifr],[nx ny inf 1]);
edia9 = ncread(ename,'ek_omb',[2 2 1 ifr],[nx ny inf 1]);

edia = edia1 + edia2 + edia3 + edia4 + edia5 + edia6 + edia7 + edia8;

sedia1 = sum(edia1.*dzhlf,3);
sedia8 = sum(edia8.*dzhlf,3);
sedia9 = sum(edia9.*dzhlf,3);

dt = t2-t1;
de = (e2-e1)/dt;

%dKe/dt = u*dia/dz

de_err = abs(edia - de);
max(de_err(:))

return




udpx = uhlf.*ncread(dname,'u_pgr',[1 1 1 1],[inf inf inf 1]);
uzxb = uhlf.*ncread(dname,'u_zxb',[1 1 1 1],[inf inf inf 1]);

p = ncread(dname,'w_prs',[1 1 1 1],[inf inf 32 1]);

zhlf = 0.5*(z1+z2);

[w,omega] = calc_w(uhlf,vhlf,zr,dz,dx,dy);
[nx,ny,nzp] = size(omega);
nz = nzp-1;
om= zeros(66,66,32);
om(2:end-1,2:end-1,:) = omega(:,:,2:end);



% d(up)
%     = u(2:end  ,2:end-1,:).*(p(3:end  ,2:end-1,:) - p(2:end-1,2:end-1,:))
%     - u(1:end-1,2:end-1,:).*(p(2:end-1,2:end-1,:) - p(1:end-2,2:end-1,:))
%   + om(2:end-1,2:end-1,0:nz).*(p(2:end-1,2:end-1,2:nz+1)-p(2:end-1,2:end-1,1:nz))
%   + om(2:end-1,2:end-1,0:nz).*(p(2:end-1,2:end-1,1:nz)-p(2:end-1,2:end-1,:nz-1))

up = uhlf(:,2:end-1,:).*(p(2:end,2:end-1,:)-p(1:end-1,2:end-1,:));
omp = zeros(nx,ny,nz+1);
omp(:,:,2:nz) = omega(:,:,2:nz).*(p(2:end-1,2:end-1,2:nz)-p(2:end-1,2:end-1,1:nz-1));

up = up.*0.5.*(dz(2:end,2:end-1,:)+dz(1:end-1,2:end-1,:));
dx = repmat(dx(2:end-1,2:end-1),[1 1 nz+1]);
omp = omp.*dx;

dup =  up(2:end,:,:) -  up(1:end-1,:,:) ...
    + omp(:,:,2:end) - omp(:,:,1:end-1);

wdzb = om.*ncread(dname,'w_dzb',[1 1 1 1],[inf inf 32 1]);
return

dzb = ncread(dname,'w_dzb',[1 1 1 1],[inf inf 32 1]);

wdia8 = sum(om.*ncread(dname,'w_dzb',[1 1 1 1],[inf inf 32 1]),3);

dia_time = ncread(dname,'ocean_time',[1],[1]);

tm1 = ncread(hnam1,'temp',[1 1 1 1],[inf inf inf 1]);
tm2 = ncread(hnam1,'temp',[1 1 1 2],[inf inf inf 1]);
tmhlf = 0.5*(tm1+tm2);

wb = sum(w.*tmhlf,3);
return

% when reading from restart files, take 2nd time index
u = ncread(hnam1,'u',[1 1 1 1],[inf inf inf 1]);
v = ncread(hnam1,'v',[1 1 1 1],[inf inf inf 1]);
zeta = ncread(hnam1,'zeta',[1 1 1],[inf inf 1]);
t1 = ncread(hnam1,'ocean_time',[1],[1]);

zw = zlevs(h,zeta,hnam1,'w');
dz = zw(:,:,2:end)-zw(:,:,1:end-1);
dzu = 0.5*(dz(2:end,:,:) + dz(1:end-1,:,:));
dzv = 0.5*(dz(:,2:end,:) + dz(:,1:end-1,:));

u1  = sum(dzu.*u,3);
v1  = sum(dzv.*v,3);
u21 = 0.5*sum(dzu.*u.*u,3);
v21 = 0.5*sum(dzv.*v.*v,3);

u = ncread(hnam2,'u',[1 1 1 2],[inf inf inf 1]);
v = ncread(hnam2,'v',[1 1 1 2],[inf inf inf 1]);
zeta = ncread(hnam2,'zeta',[1 1 2],[inf inf 1]);
t2 = ncread(hnam2,'ocean_time',[2],[1]);

zw = zlevs(h,zeta,hnam2,'w');
dz = zw(:,:,2:end)-zw(:,:,1:end-1);
dzu = 0.5*(dz(2:end,:,:) + dz(1:end-1,:,:));
dzv = 0.5*(dz(:,2:end,:) + dz(:,1:end-1,:));

u2 = sum(dzu.*u,3);
v2 = sum(dzv.*v,3);
u22 = 0.5*sum(dzu.*u.*u,3);
v22 = 0.5*sum(dzv.*v.*v,3);

dv2 = v22-v21;
du2 = u22-u21;
du = u2-u1;
dt = t2-t1;

udia = udia1+udia2+udia3+udia4+udia5+udia6+udia7;
vdia = vdia1+vdia2+vdia3+vdia4+vdia5+vdia6+vdia7;

resu = udia-du2/dt;
resv = vdia-dv2/dt;

edia = u2rho(udia) + v2rho(vdia);
de2  = u2rho(du2) + v2rho(dv2);
[max(du2(:)) max(dv2(:)) max(de2(:)) ]
[du2(14,16) dv2(14,16)]

rese =  edia-de2/dt;


resu(1,:) = 0; resu(end,:) = 0;
resu(:,1) = 0; resu(:,end) = 0;
resv(1,:) = 0; resv(end,:) = 0;
resv(:,1) = 0; resv(:,end) = 0;
rese(1,:) = 0; rese(end,:) = 0;
rese(:,1) = 0; rese(:,end) = 0;

imagesc(rese');axis xy;colorbar
cmocean('balance');










