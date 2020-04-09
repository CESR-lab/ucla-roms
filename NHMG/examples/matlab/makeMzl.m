function M = makeMzl(dx,dy,dz,nx,ny,nz);

  % zlevel operator for constant grid spacing
  
  ndu = (nx+1)*ny*nz;
  ndv = nx*(ny+1)*nz;
  ndw = nx*ny*(nz+1);

  cuc = dx*dy*dz*ones(nx+1,ny,nz);
  cvc = dx*dy*dz*ones(nx,ny+1,nz);
  cwc = dx*dy*dz*ones(nx,ny,nz+1);
  cwc(:,:,nz+1) = 0.5*cwc(:,:,nz+1);
  
  % Put coefficients in matrix
  nd = ndu + ndv + ndw;
  iM = zeros(nd,1);
  jM = zeros(nd,1);
  cM = zeros(nd,1); 
   
  v = 1;
   
  %% Coefficients for u
  [i,j,k] = meshgrid([1:nx+1],[1:ny],[1:nz]);
  i = reshape(permute(i,[2 1 3]),ndu,1);
  j = reshape(permute(j,[2 1 3]),ndu,1);
  k = reshape(permute(k,[2 1 3]),ndu,1);   
  udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
  
  [iM,jM,cM,v] = addcoef(iM,jM,cM,v,cuc,udx,udx);


  %% Coefficients for vf
  [i,j,k] = meshgrid([1:nx],[1:ny+1],[1:nz]);
  i = reshape(permute(i,[2 1 3]),ndv,1);
  j = reshape(permute(j,[2 1 3]),ndv,1);
  k = reshape(permute(k,[2 1 3]),ndv,1);      
  vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
  
  [iM,jM,cM,v] = addcoef(iM,jM,cM,v,cvc,vdx,vdx);

   
  %% Coefficients for wf
  [i,j,k] = meshgrid([1:nx],[1:ny],[1:nz+1]);
  i = reshape(permute(i,[2 1 3]),ndw,1);
  j = reshape(permute(j,[2 1 3]),ndw,1);
  k = reshape(permute(k,[2 1 3]),ndw,1);  
  wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;

  [iM,jM,cM,v] = addcoef(iM,jM,cM,v,cwc,wdx,wdx);
  
  M = sparse(iM,jM,cM,nd,nd,nd);

end

