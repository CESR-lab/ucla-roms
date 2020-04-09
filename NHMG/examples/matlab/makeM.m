function M = makeM(Au,Av,Aw,dxu,dyv,dzw)

[nxp,ny,nz] = size(Au);
nx = nxp-1;
ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);

if 0
    
M = sparse([],[],[],ndu+ndv+ndw,ndu+ndv+ndw,ndu+ndv+ndw);
for i = 1:nx+1
  for j = 1:ny  
    for k = 1:nz
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;        
      M(udx+1,udx+1) = Au(i,j,k)*dxu(i,j);
    end
  end
end

for i = 1:nx
  for j = 1:ny+1  
    for k = 1:nz
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      M(vdx+1,vdx+1) = Av(i,j,k)*dyv(i,j);
    end
  end
end

for i = 1:nx
  for j = 1:ny  
    for k = 1:nz+1
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i-1 + ndu +ndv;
      M(wdx+1,wdx+1) = Aw(i,j)*dzw(i,j,k);
    end
  end
end
else
  cuc = Au.*repmat(dxu,[1 1 nz]);
  cvc = Av.*repmat(dyv,[1 1 nz]);
  cwc = repmat(Aw,[1 1 nz+1]).*dzw;
  
  %% Put coefficients in matrix
  nd = ndu + ndv + ndw;
  iM = zeros(nd,1);
  jM = zeros(nd,1);
  cM = zeros(nd,1); 
   
  idx = 1;
   
  %% Coefficients for u
  [i,j,k] = meshgrid([1:nx+1],[1:ny],[1:nz]);
  i = reshape(permute(i,[2 1 3]),ndu,1);
  j = reshape(permute(j,[2 1 3]),ndu,1);
  k = reshape(permute(k,[2 1 3]),ndu,1);   
  udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
   
  iM(idx:idx+ndu-1) = udx;
  jM(idx:idx+ndu-1) = udx;
  cM(idx:idx+ndu-1) = reshape(cuc,ndu,1);
  idx = idx + ndu;

  %% Coefficients for vf
  [i,j,k] = meshgrid([1:nx],[1:ny+1],[1:nz]);
  i = reshape(permute(i,[2 1 3]),ndv,1);
  j = reshape(permute(j,[2 1 3]),ndv,1);
  k = reshape(permute(k,[2 1 3]),ndv,1);      
  vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;

  iM(idx:idx+ndv-1) = vdx;
  jM(idx:idx+ndv-1) = vdx;
  cM(idx:idx+ndv-1) = reshape(cvc,ndv,1);
  idx = idx + ndv;  
   
  %% Coefficients for wf
  [i,j,k] = meshgrid([1:nx],[1:ny],[1:nz+1]);
  i = reshape(permute(i,[2 1 3]),ndw,1);
  j = reshape(permute(j,[2 1 3]),ndw,1);
  k = reshape(permute(k,[2 1 3]),ndw,1);  
  wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;

  iM(idx:idx+ndw-1) = wdx;
  jM(idx:idx+ndw-1) = wdx;
  cM(idx:idx+ndw-1) = reshape(cwc,ndw,1);
  idx = idx + ndw;  
  
  M = sparse(iM,jM,cM,nd,nd,nd);
  clear 'iM' 'jM' 'cM'

end

