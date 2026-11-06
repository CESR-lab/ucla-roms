% This code checks the online momentum diagnostics. 
% There are currently 2 versions of the diagnostics;
% normalized and unnormalized, where the diagnostic
% terms include a multiplication with dz. 
% The normalized version should do better in capturing
% the free surface term, u dz/dt /dz, which could be
% significant in shallow tidal areas. 
% However, the non-normalized version does better 
% computing the wb term. 


gname = 'input_data/sample_grd_riv.nc';
fname = 'diags_his.20121209133435.nc';
dname = 'diags_dia.20121209143435.nc';

ifr = 1;

% diagnostics can include at multiplication by dz or not (diag_norm = 0/1)
if findstr(ncreadatt(dname,'/','normalization'),'include');
  diag_norm = logical(0);
else
  diag_norm = logical(1);
end

[nxu,nyp,nz] = size(ncread(fname,'u',[1 1 1 ifr],[inf inf inf 1]));
nx = nxu-1;
ny = nyp-2;
nyv= nyp-1;

h    = ncread(gname,'h');
msk  = ncread(gname,'mask_rho');

msku = msk(2:end,2:end-1).*msk(1:end-1,2:end-1);
mskv = msk(2:end-1,2:end).*msk(2:end-1,1:end-1);
msku3d = repmat(msku,[1 1 nz]);
mskv3d = repmat(mskv,[1 1 nz]);

tim = ncread(fname,'ocean_time',[1],[2]);
delt = tim(2:end)-tim(1:end-1);
delt = sum(delt);

z1= ncread(fname,'zeta',[1 1 ifr],[inf inf 1]);
zw1= zlevs(h,z1,fname,'w');
dz1 = zw1(:,:,2:end)- zw1(:,:,1:end-1); 
dzu1 = 0.5*(dz1(2:end,2:end-1,:) + dz1(1:end-1,2:end-1,:));
dzv1 = 0.5*(dz1(2:end-1,2:end,:) + dz1(2:end-1,1:end-1,:));
t1 = ncread(fname,'ocean_time',[ifr],[1]);
u1 = ncread(fname,'u',[1 2 1 ifr],[nxu ny inf 1]);
v1 = ncread(fname,'v',[2 1 1 ifr],[nx nyv inf 1]);

z2= ncread(fname,'zeta',[1 1 ifr+1],[inf inf 1]);
zw2= zlevs(h,z2,fname,'w');
dz2 = zw2(:,:,2:end)- zw2(:,:,1:end-1); 
dzu2 = 0.5*(dz2(2:end,2:end-1,:) + dz2(1:end-1,2:end-1,:));
dzv2 = 0.5*(dz2(2:end-1,2:end,:) + dz2(2:end-1,1:end-1,:));
t2 = ncread(fname,'ocean_time',[ifr+1],[1]);
u2 = ncread(fname,'u',[1 2 1 ifr+1],[nxu ny inf 1]);
v2 = ncread(fname,'v',[2 1 1 ifr+1],[nx nyv inf 1]);

% The boundary conditions are not included in the diagnostics
u1(1,:,:) = 0; u1(end,:,:) = 0;
u2(1,:,:) = 0; u2(end,:,:) = 0;
v1(:,1,:) = 0; v1(:,end,:) = 0;
v2(:,1,:) = 0; v2(:,end,:) = 0;

if diag_norm
  udz1 = u1;
  udz2 = u2;
else
  udz1 = dzu1.*u1;
  udz2 = dzu2.*u2;
end

deludz = udz2-udz1;

 urhs1 = ncread(dname,'u_pgr',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs2 = ncread(dname,'u_cor',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs3 = ncread(dname,'u_adv',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs4 = ncread(dname,'u_dis',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs5 = ncread(dname,'u_hmx',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs6 = ncread(dname,'u_vmx',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs7 = ncread(dname,'u_cpl',[1 2 1 1],[nxu ny inf 1]).*msku3d;
 urhs8 = ncread(dname,'u_dzt',[1 2 1 1],[nxu ny inf 1]).*msku3d;

 if diag_norm
   urhs = urhs1 + urhs2 + urhs3 + urhs4 + urhs5 + urhs6 + urhs7 + urhs8;
 else
   urhs = urhs1 + urhs2 + urhs3 + urhs4 + urhs5 + urhs6 + urhs7;
 end


 urhs = delt*urhs;

 u_err = (deludz - urhs);
 max(abs(u_err(:)))


  imagesc(u_err(:,:,18)');axis xy;colorbar
