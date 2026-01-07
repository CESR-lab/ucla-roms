gname = 'input_data/sample_grd_riv.nc';
dname = 'diags_dia.20121209143435.nc';
hname = 'diags_his.20121209133435.nc';
ename = 'diags_ekb.20121209143435.nc';

h = ncread(gname,'h');
h = h(2:end-1,2:end-1);
[nx,ny] = size(h);
h = ncread(gname,'h');
nxu = nx+1;
nyv = ny+1;
dx = 1./ncread(gname,'pm');
dy = 1./ncread(gname,'pn');

for ifr= 1:1
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
zr1 = zlevs(h,z1,hname,'r');
dz1 = zw1(:,:,2:end)-zw1(:,:,1:end-1);

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


edia1 = ncread(ename,'ek_pgr',[2 2 ifr],[nx ny 1]);
edia2 = ncread(ename,'ek_cor',[2 2 ifr],[nx ny 1]);
edia3 = ncread(ename,'ek_adv',[2 2 ifr],[nx ny 1]);
edia4 = ncread(ename,'ek_dis',[2 2 ifr],[nx ny 1]);
edia5 = ncread(ename,'ek_hmx',[2 2 ifr],[nx ny 1]);
edia6 = ncread(ename,'ek_vmx',[2 2 ifr],[nx ny 1]);
edia7 = ncread(ename,'ek_cpl',[2 2 ifr],[nx ny 1]);
edia8 = ncread(ename,'ek_dzt',[2 2 ifr],[nx ny 1]);
edia9 = ncread(ename,'ek_omb',[2 2 ifr],[nx ny 1]);

edia = edia1 + edia2 + edia3 + edia4 + edia5 + edia6 + edia7 + edia8;

e1 = sum(e1.*dz1(2:end-1,2:end-1,:),3);
e2 = sum(e2.*dz2(2:end-1,2:end-1,:),3);

dt = t2-t1;
de = (e2-e1)/dt; %./dzhlf;

%dKe/dt = u*dia/dz

de_err = abs(edia - de);
max(de_err(:))


whlf = 0.5*(ncread(hname,'w',[2 2 1 ifr],[nx ny inf 1])+ncread(hname,'w',[2 2 1 ifr+1],[nx ny inf 1]));
bhlf = 9.81*0.5*(ncread(hname,'rho',[2 2 1 ifr],[nx ny inf 1])+ncread(hname,'rho',[2 2 1 ifr+1],[nx ny inf 1]))/1e3;
dzhlf = 0.5*(dz1+dz2);
dzhlf = dzhlf(2:end-1,2:end-1,:);

wb = whlf.*bhlf;
wb = sum(wb.*dzhlf,3);

swb(ifr) = sum(wb(:));

sadv(ifr) = sum(edia3(:));
scpl(ifr) = sum(edia7(:));
sdzt(ifr) = sum(edia8(:));
spgr(ifr) = sum(edia1(:));
somb(ifr) = sum(edia9(:));

end
