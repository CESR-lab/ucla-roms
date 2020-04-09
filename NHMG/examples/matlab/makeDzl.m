function D = makeDzl(dx,dy,dz,mask)

   [nx,ny,nz] = size(mask);
   nd = nx*ny*nz;
   ndu = (nx+1)*ny*nz;
   ndv = nx*(ny+1)*nz;
   ndw = nx*ny*(nz+1);

   cuc = zeros(nx,ny,nz);
   cue = zeros(nx,ny,nz);
   cvc = zeros(nx,ny,nz);
   cvn = zeros(nx,ny,nz);
   cwc = zeros(nx,ny,nz);
   cwu = zeros(nx,ny,nz);  

   cuc(2:nx  ,:,:) = -dy*dz.*mask(2:nx,:,:).*mask(1:nx-1,:,:);  
   cue(1:nx-1,:,:) =  dy*dz.*mask(2:nx,:,:).*mask(1:nx-1,:,:);
   cvc(:,2:ny  ,:) = -dx*dz.*mask(:,2:ny,:).*mask(:,1:ny-1,:);
   cvn(:,1:ny-1,:) =  dx*dz.*mask(:,2:ny,:).*mask(:,1:ny-1,:);
   cwc(:,:,2:nz  ) = -dx*dy.*mask(:,:,2:nz).*mask(:,:,1:nz-1);
   cwu(:,:,1:nz-1) =  dx*dy.*mask(:,:,2:nz).*mask(:,:,1:nz-1); 
  % max(abs(cvc))
   
   cwu(:,:,nz) = dx*dy.*mask(:,:,nz); 
   
   cuc(mask<1) = 1;
   

    
   %% Put coefficients in matrix
       
   iD = zeros(6*nd,1);
   jD = zeros(6*nd,1);
   cD = zeros(6*nd,1); 
   v = 1;
   
   [i,j,k] = meshgrid([1:nx],[1:ny],[1:nz]);
   i = reshape(permute(i,[2 1 3]),nd,1);
   j = reshape(permute(j,[2 1 3]),nd,1);
   k = reshape(permute(k,[2 1 3]),nd,1);   
   idx = (k-1)*ny*nx + (j-1)*nx + i;
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i; 
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
   wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
     
   is = 1;
   js = nx;
   ks = nx*ny;
   
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cuc,idx,udx);
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cue,idx,udx+is);   
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cvc,idx,vdx);
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cvn,idx,vdx+js);   
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cwc,idx,wdx);
   [iD,jD,cD,v] = addcoef(iD,jD,cD,v,cwu,idx,wdx+ks);      

   iD(cD==0) = [];
   jD(cD==0) = [];
   cD(cD==0) = [];
   nnz = length(cD);     
   
   D = sparse(iD,jD,cD,nd,ndu+ndv+ndw,nnz);

end


