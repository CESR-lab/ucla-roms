function D = makeD(Au,Av,Aw,rname)

[nxp,ny,nz] = size(Au);
nx = nxp-1;
nd = nx*ny*nz;
ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);

if strcmp(rname,'Cbc')
  ndb = nx*ny;
else
  ndb=0;
end


is = 1;
js = nx;
ks = nx*ny;

if 0
D = sparse(nd+nx*ny,ndu+ndv+ndw,6*nd);
for i = 1:nx
  for j = 1:ny  
    for k = 1:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i; 
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
      if i>1
        D(idx,udx   ) = -Au(i  ,j,k);
      end
      if i<nx
        D(idx,udx+is) =  Au(i+1,j,k);
      end
      if j>1
        D(idx,vdx   ) = -Av(i,j  ,k);
      end
      if j<ny
        D(idx,vdx+js) =  Av(i,j+1,k);
      end
      if k>1
        D(idx,wdx   ) = -Aw(i,j);
      end
      D(idx,wdx+ks) =    Aw(i,j);

    end
  end
end

 disp('extended D')
 for i = 1:nx
  for j = 1:ny
    k=1;
    idx = (j-1)*nx + i +nd;
    wdx = (k-1)*ny*nx     + (j-1)*nx    + i + ndu +ndv;
    D(idx,wdx) =    Aw(i,j);    
  end
 end

else
   disp('Fast makeD')
   tic
   cuc = zeros(nx,ny,nz);
   cue = zeros(nx,ny,nz);
   cvc = zeros(nx,ny,nz);
   cvn = zeros(nx,ny,nz);
   cwc = zeros(nx,ny,nz);
   cwu = zeros(nx,ny,nz);  
   
   cuc(2:nx  ,:,:) = -Au(2:nx,:,:);  
   cue(1:nx-1,:,:) =  Au(2:nx,:,:);
   cvc(:,2:ny  ,:) = -Av(:,2:ny,:);
   cvn(:,1:ny-1,:) =  Av(:,2:ny,:);
   cwc(:,:,2:nz  ) = -repmat(Aw,[1 1 nz-1]);
   cwu(:,:,1:nz  ) =  repmat(Aw,[1 1 nz]); 
    
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
   
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = udx;
   cD(v:v+nd-1) = reshape(cuc,nd,1);
   v = v + nd;
   
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = udx +is;
   cD(v:v+nd-1) = reshape(cue,nd,1);
   v = v + nd;
   
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = vdx;
   cD(v:v+nd-1) = reshape(cvc,nd,1);
   v = v + nd;
   
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = vdx +js;
   cD(v:v+nd-1) = reshape(cvn,nd,1);
   v = v + nd;
    
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = wdx;
   cD(v:v+nd-1) = reshape(cwc,nd,1);
   v = v + nd;
   
   iD(v:v+nd-1) = idx;
   jD(v:v+nd-1) = wdx +ks;
   cD(v:v+nd-1) = reshape(cwu,nd,1);
   v = v + nd;
    

   % extended D
   if ndb>0
    disp('extended D')
    [i,j] = meshgrid([1:nx],[1:ny]);
    i = reshape(permute(i,[2 1]),ndb,1);
    j = reshape(permute(j,[2 1]),ndb,1);
    idx = (j-1)*nx + i + nd;
    wdx = (j-1)*nx     + i + ndu +ndv;
   
    iD(v:v+ndb-1) = idx;
    jD(v:v+ndb-1) = wdx;
    cD(v:v+ndb-1) = 1;
    v = v + ndb;
   end
   
   iD(cD==0) = [];
   jD(cD==0) = [];
   cD(cD==0) = [];
   nnz = length(cD);     
   
   D = sparse(iD,jD,cD,nd+ndb,ndu+ndv+ndw,nnz);
   clear 'iD' 'jD' 'cD'
   toc
end


