function Tc = makeTc(zxdy,zydx,dxu,dyv,Aw)

[nx,ny,nz] = size(zxdy);
ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);

% Flux representation: (uf,vf,wf)

% Fom Flux to real
% ur = uf
% vr = vf
% wr = uf*zx + vf*zy + wf

% Real to Flux
% uf = ur
% vf = vr
% wf = wr - ur*zx - vr*zy


% Ek in terms of uf,wf
%
%  Ek = 0.5*uf^2  + 0.5*vf^2 + 0.5*(uf*zx + vf*zy + wf)^2


% find (um,vm,wm) such that Ek = umuf + vmvf + wmwf

%  um = Ek_uf = uf*(1+zx^2) + vf*zx*zy + wf*zx 
%  vm = Ek_vf = vf*(1+zy^2) + uf*zx*zy + wf*zy
%  wm = Ek_wf = uf*zx +vf*zy   + wf

%  Invert this to find:
%
%  uf =  um - wm*zx
%  vf =  vm - wm*zy
%%  wf = -um*z_x - vm*zy +wm*(1+zy^2+ zx^2)
%zxdy = 0*zxdy;
%zydx = 0*zydx;
%%
if 0
tic
disp('slow makeTc')
Tc = sparse(ndu+ndv+ndw,ndu+ndv+ndw,ndu+ndv+9*ndw);
is = 1;
js = nx;
ks = nx*ny;
for i = 1:nx+1
  for j = 1:ny  
    for k = 1:nz
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;  
      %  uf =  uc
      Tc(udx+1,udx+1) = 1; 
    end
  end 
end

jsu = nx+1;
for i = 1:nx
  for j = 1:ny+1
    for k = 1:nz 
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      %  vf =  vc 
      Tc(vdx+1,vdx+1) = 1;             
    end
  end 
end

ksu = ny*(nx+1);
ksv = (ny+1)*nx;
% wf = -zx*uc+-zy*vc + *wc
for i = 1:nx
  for j = 1:ny    
    for k = 1:nz+1
        
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;  
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i-1 + ndu +ndv;      
      
      Tc(wdx+1,wdx+1)= 1;
      if k<nz+1&k>1

        Tc(wdx+1,udx+1   -ksu)=  -0.25*zxdy(i,j,k-1)*dxu(i  ,j)/Aw(i,j);
        Tc(wdx+1,udx+1       )=  -0.25*zxdy(i,j,k  )*dxu(i  ,j)/Aw(i,j);
        Tc(wdx+1,udx+1+is-ksu)=  -0.25*zxdy(i,j,k-1)*dxu(i+1,j)/Aw(i,j);
        Tc(wdx+1,udx+1+is    )=  -0.25*zxdy(i,j,k  )*dxu(i+1,j)/Aw(i,j);
        
        Tc(wdx+1,vdx+1   -ksv)=  -0.25*zydx(i,j,k-1)*dyv(i,j  )/Aw(i,j);
        Tc(wdx+1,vdx+1       )=  -0.25*zydx(i,j,k  )*dyv(i,j  )/Aw(i,j);
        Tc(wdx+1,vdx+1+js-ksv)=  -0.25*zydx(i,j,k-1)*dyv(i,j+1)/Aw(i,j);
        Tc(wdx+1,vdx+1+js    )=  -0.25*zydx(i,j,k  )*dyv(i,j+1)/Aw(i,j);
        
      elseif k==nz+1  
  
        Tc(wdx+1,udx+1   -ksu)=  -0.5*zxdy(i,j,k-1)*dxu(i  ,j)/Aw(i,j);
        Tc(wdx+1,udx+1+is-ksu)=  -0.5*zxdy(i,j,k-1)*dxu(i+1,j)/Aw(i,j);
        
        Tc(wdx+1,vdx+1   -ksv)=  -0.5*zydx(i,j,k-1)*dyv(i,j  )/Aw(i,j);
        Tc(wdx+1,vdx+1+js-ksv)=  -0.5*zydx(i,j,k-1)*dyv(i,j+1)/Aw(i,j);
        
      elseif k==1

 %       Tc(wdx+1,udx+1   +ksu)=  0.25*zxdy(i,j,k+1)*dxu(i  ,j)/Aw(i,j);
 %       Tc(wdx+1,udx+1+is+ksu)=  0.25*zxdy(i,j,k+1)*dxu(i+1,j)/Aw(i,j);
        
        Tc(wdx+1,udx+1       )=  -0.5*zxdy(i,j,k  )*dxu(i  ,j)/Aw(i,j);
        Tc(wdx+1,udx+1+is    )=  -0.5*zxdy(i,j,k  )*dxu(i+1,j)/Aw(i,j);
        
        Tc(wdx+1,vdx+1       )=  -0.5*zydx(i,j,k  )*dyv(i,j  )/Aw(i,j);
        Tc(wdx+1,vdx+1+js    )=  -0.5*zydx(i,j,k  )*dyv(i,j+1)/Aw(i,j);
      else
         error('kan niet')
          
      end

    end      
  end
end
toc
else
   disp('fast makeTc')
   tic

   cuucc = ones(nx+1,ny,nz);
   cvvcc = ones(nx,ny+1,nz);  
   cwwcc = ones(nx,ny,nz+1);
   
   cwucc = zeros(nx,ny,nz+1);
   cwuec = zeros(nx,ny,nz+1);
   cwucd = zeros(nx,ny,nz+1);
   cwued = zeros(nx,ny,nz+1);
   cwvcc = zeros(nx,ny,nz+1);
   cwvnc = zeros(nx,ny,nz+1);
   cwvcd = zeros(nx,ny,nz+1);
   cwvnd = zeros(nx,ny,nz+1);   

   cwucc(:,:,1:nz  ) = -0.25*zxdy.*repmat(dxu(1:nx  ,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwuec(:,:,1:nz  ) = -0.25*zxdy.*repmat(dxu(2:nx+1,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwucd(:,:,2:nz+1) = -0.25*zxdy.*repmat(dxu(1:nx  ,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwued(:,:,2:nz+1) = -0.25*zxdy.*repmat(dxu(2:nx+1,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   
   cwvcc(:,:,1:nz  ) = -0.25*zydx.*repmat(dyv(:,1:ny  ),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvnc(:,:,1:nz  ) = -0.25*zydx.*repmat(dyv(:,2:ny+1),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvcd(:,:,2:nz+1) = -0.25*zydx.*repmat(dyv(:,1:ny  ),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvnd(:,:,2:nz+1) = -0.25*zydx.*repmat(dyv(:,2:ny+1),[1 1 nz])./repmat(Aw,[1 1 nz]); 
   
   if 0
     cwwcc(:,:,1) = 0;
     cwucc(:,:,1) = 0;cwuec(:,:,1) = 0;
     cwvcc(:,:,1) = 0;cwvnc(:,:,1) = 0;
   else
     cwucc(:,:,1) = 2*cwucc(:,:,1);
     cwuec(:,:,1) = 2*cwuec(:,:,1);
     cwvcc(:,:,1) = 2*cwvcc(:,:,1);
     cwvnc(:,:,1) = 2*cwvnc(:,:,1);
       
   end
   
   cwucd(:,:,nz+1) = 2*cwucd(:,:,nz+1);
   cwued(:,:,nz+1) = 2*cwued(:,:,nz+1);
   cwvcd(:,:,nz+1) = 2*cwvcd(:,:,nz+1);
   cwvnd(:,:,nz+1) = 2*cwvnd(:,:,nz+1);    

   %% Put coefficients in matrix
   nnz = ndu + ndu + 9*ndw;
   iT = zeros(nnz,1);
   jT = zeros(nnz,1);
   cT = zeros(nnz,1); 
   
   is = 1; js=nx;ks = nx*ny;
   idx = 1;
   
   %% Coefficients for uf
   [i,j,k] = meshgrid([1:nx+1],[1:ny],[1:nz]);
   i = reshape(permute(i,[2 1 3]),ndu,1);
   j = reshape(permute(j,[2 1 3]),ndu,1);
   k = reshape(permute(k,[2 1 3]),ndu,1);   
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
   
   iT(idx:idx+ndu-1) = udx;
   jT(idx:idx+ndu-1) = udx;
   cT(idx:idx+ndu-1) = reshape(cuucc,ndu,1);
   idx = idx + ndu;
   
   %% Coefficients for vf
   [i,j,k] = meshgrid([1:nx],[1:ny+1],[1:nz]);
   i = reshape(permute(i,[2 1 3]),ndv,1);
   j = reshape(permute(j,[2 1 3]),ndv,1);
   k = reshape(permute(k,[2 1 3]),ndv,1);      
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;

   iT(idx:idx+ndv-1) = vdx;
   jT(idx:idx+ndv-1) = vdx;
   cT(idx:idx+ndv-1) = reshape(cvvcc,ndv,1);
   idx = idx + ndv;  
   
   %% Coefficients for wf
   [i,j,k] = meshgrid([1:nx],[1:ny],[1:nz+1]);
   i = reshape(permute(i,[2 1 3]),ndw,1);
   j = reshape(permute(j,[2 1 3]),ndw,1);
   k = reshape(permute(k,[2 1 3]),ndw,1);  
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
   wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
   is = 1;
   ksu = ny*(nx+1);
   ksv = (ny+1)*nx;

   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = wdx;
   cT(idx:idx+ndw-1) = reshape(cwwcc,ndw,1);
   idx = idx + ndw;  
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = udx;
   cT(idx:idx+ndw-1) = reshape(cwucc,ndw,1);
   idx = idx + ndw;  
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = udx+is;
   cT(idx:idx+ndw-1) = reshape(cwuec,ndw,1);
   idx = idx + ndw;  
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = udx-ksu;
   cT(idx:idx+ndw-1) = reshape(cwucd,ndw,1);
   idx = idx + ndw; 
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = udx+is-ksu;
   cT(idx:idx+ndw-1) = reshape(cwued,ndw,1);
   idx = idx + ndw; 
    
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = vdx;
   cT(idx:idx+ndw-1) = reshape(cwvcc,ndw,1);
   idx = idx + ndw;  
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = vdx+js;
   cT(idx:idx+ndw-1) = reshape(cwvnc,ndw,1);
   idx = idx + ndw;  
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = vdx-ksv;
   cT(idx:idx+ndw-1) = reshape(cwvcd,ndw,1);
   idx = idx + ndw; 
   
   iT(idx:idx+ndw-1) = wdx;
   jT(idx:idx+ndw-1) = vdx+js-ksv;
   cT(idx:idx+ndw-1) = reshape(cwvnd,ndw,1);
   idx = idx + ndw;     
   
   iT(cT==0) = [];
   jT(cT==0) = [];
   cT(cT==0) = [];
   nnz = length(cT);
   
   Tc = sparse(iT,jT,cT,ndu+ndv+ndw,ndu+ndv+ndw,nnz);
   clear 'iT' 'jT' 'cT'
toc
end