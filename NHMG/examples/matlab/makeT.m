function T = makeT(zx,zy,zxdy,zydx,Au,Av,Aw,dxu,dyv,dz,dzw,rname)

%% Flux representation: (uf,vf,wf)

% Fom Flux to real
% ur = uf
% vr = vf
% wr = uf*zx + vf*zy + wf

%% Real to Flux
% uf = ur
% vf = vr
% wf = wr - ur*zx - vr*zy

%% Real to momentum
% um = ur + wr*zx
% vm = vr + wr*zy
% wm = wr 

%% Momentum to real
% ur = um = wm*zx
% vr = vm - wm*zy
% wr = wm

%% Ek in terms of uf,wf
%
%  Ek = 0.5*uf^2  + 0.5*vf^2 + 0.5*(uf*zx + vf*zy + wf)^2


%% find (um,vm,wm) such that Ek = umuf + vmvf + wmwf
%  um = Ek_uf = uf*(1+zx^2) + vf*zx*zy + wf*zx 
%  vm = Ek_vf = vf*(1+zy^2) + uf*zx*zy + wf*zy
%  wm = Ek_wf = uf*zx +vf*zy   + wf

%%  Invert this to find:
%
%  uf =  um - wm*zx
%  vf =  vm - wm*zy
%  wf = -um*z_x - vm*zy +wm*(1+zy^2+ zx^2)


%%wm = um*zx - wm*zx^2 + vm*zy - wm*zy^2 -um*zx - vm*zy + wm*(1+zx^2+zy^2)
%%um = um*(1+zx^2) - wm*zx*(1+zx^2) + vm*zx*zy - wm*zx*zy^2 - um*zx^2 -vm*zy*zx + wm*zx*(1+zx^2 + zy^2)
%%
%% Unbelievable, but true :)

[nx,ny,nz] = size(zx);
ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);

alpha = (1 + zx.^2 + zy.^2); %% defined at rho-points
zxzy8dzcwi = 0.125*zx(:,:,1).*zy(:,:,1).*dz(:,:,1)./alpha(:,:,1);

if 0
    tic
T = sparse([],[],[],ndu+ndv+ndw,ndu+ndv+ndw,5*ndu+5*ndv+9*ndw);
is = 1;
js = nx;
ks = nx*ny;
for i = 1:nx
  for j = 1:ny  
    for k = 1:nz
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;  
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i-1 + ndu +ndv;

      %%  uf =  um - wm*zx
      T(udx+1,udx+1) = 1; 
      if i>1&&i<nx+1 %% The slope is zero at the side boundaries   
        
        T(udx+1,wdx+1-is+ks) = -0.25*zxdy(i-1,j,k)*dzw(i-1,j,k+1)/Au(i,j,k);
        T(udx+1,wdx+1   +ks) = -0.25*zxdy(i  ,j,k)*dzw(i  ,j,k+1)/Au(i,j,k);
        
        if k>1
          T(udx+1,wdx+1-is   ) = -0.25*zxdy(i-1,j,k)*dzw(i-1,j,k  )/Au(i,j,k);
          T(udx+1,wdx+1      ) = -0.25*zxdy(i  ,j,k)*dzw(i  ,j,k  )/Au(i,j,k);
        else
          %% W = -u*zx - v*zy - (1+zx^2+zy^2)w = 0 at z=0
          % U(i,j,1) = u(i,j,1) - 0.25*w(i-1,j,2)*zx(i-1,j,1) - 0.25*w(i,j,2)*zx(i,j,1)
          %                     - 0.5*w_u(i,j
          %  w_u =    zx_u*u/alpha_u +
          %     + 0.25*( (zy(i-1,j,1)*(v(i-1,j,1)+ v(i-1,j+1,1))/cw(i-1,j)
          %            + (zy(i  ,j,1)*(v(i  ,j,1)+ v(i  ,j+1,1))/cw(i  ,j) )
          if 1
            T(udx+1,udx+1) = 1 - 0.25*(zx(i-1,j,1)^2*cwi(i-1,j)+zx(i,j,1)^2*cwi(i,j)); 
            T(udx+1,vdx+1-is   ) = -zxzy8dzcwi(i-1,j)*dyv(i-1,j  )/Au(i,j,k);
            T(udx+1,vdx+1-is+js) = -zxzy8dzcwi(i-1,j)*dyv(i-1,j+1)/Au(i,j,k);
            T(udx+1,vdx+1      ) = -zxzy8dzcwi(i  ,j)*dyv(i  ,j  )/Au(i,j,k);
            T(udx+1,vdx+1   +js) = -zxzy8dzcwi(i  ,j)*dyv(i  ,j+1)/Au(i,j,k);
          else
            T(udx+1,wdx+1-is+ks) = 2*T(udx+1,wdx+1-is+ks);
            T(udx+1,wdx+1   +ks) = 2*T(udx+1,wdx+1   +ks);
          end
        end
        if k==nz
          T(udx+1,wdx+1-is+ks) = 2*T(udx+1,wdx+1-is+ks);
          T(udx+1,wdx+1   +ks) = 2*T(udx+1,wdx+1   +ks);
        end

      end
      
    end
  end 
end

jsu = nx+1;
for i = 1:nx
  for j = 1:ny
    for k = 1:nz
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;  
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i-1 + ndu +ndv;

        
      T(vdx+1,vdx+1) = 1;         %%  vf =  vm - wm*zy
      if j>1&&j<ny+1 %% The slope is zero at the side boundaries
        T(vdx+1,wdx+1-js+ks) = -0.25*zydx(i,j-1,k)*dzw(i,j-1,k+1)/Av(i,j,k);
        T(vdx+1,wdx+1   +ks) = -0.25*zydx(i,j  ,k)*dzw(i,j  ,k+1)/Av(i,j,k);
        if k>1
          T(vdx+1,wdx+1-js   ) = -0.25*zydx(i,j-1,k)*dzw(i,j-1,k  )/Av(i,j,k);
          T(vdx+1,wdx+1      ) = -0.25*zydx(i,j  ,k)*dzw(i,j  ,k  )/Av(i,j,k);
        else
          %% W = -u*zx - v*zy - (1+zx^2+zy^2)w = 0 at z=0
          % V(i,j,1) = v(i,j,1) - 0.25*w(i,j-1,2)*zy(i,j-1,1) - 0.25*w(i,j,2)*zy(i,j,1)
          %                     - 0.5*w_v(i,j
          %  w_v =    zy_v*v/cw +
          %     + 0.25*( (zy(i-1,j,1)*(v(i-1,j,1)+ v(i-1,j+1,1))/cw(i-1,j)
          %            + (zy(i  ,j,1)*(v(i  ,j,1)+ v(i  ,j+1,1))/cw(i  ,j) )
          if 1
            T(vdx+1,vdx+1) = 1 - 0.25*(zy(i,j-1,1)^2*cwi(i,j-1)+zy(i,j,1)^2*cwi(i,j)); 
            T(vdx+1,udx+1-jsu   ) = -0.125*zxzydzcwi(i,j-1)*dxu(i  ,j-1)/Av(i,j,k);
            T(vdx+1,udx+1-jsu+is) = -0.125*zxzydzcwi(i,j-1)*dxu(i+1,j-1)/Av(i,j,k);
            T(vdx+1,udx+1       ) = -0.125*zxzydzcwi(i  ,j)*dxu(i  ,j  )/Av(i,j,k);
            T(vdx+1,udx+1    +is) = -0.125*zxzydzcwi(i  ,j)*dxu(i+1,j  )/Av(i,j,k);
          else
            T(vdx+1,wdx+1-js+ks) = 2*T(vdx+1,wdx+1-js+ks);
            T(vdx+1,wdx+1   +ks) = 2*T(vdx+1,wdx+1   +ks);
          end
        end
        if k==nz
          T(vdx+1,wdx+1-js+ks) = 2*T(vdx+1,wdx+1-js+ks);
          T(vdx+1,wdx+1   +ks) = 2*T(vdx+1,wdx+1   +ks);
        end
      end
      
    end
  end 
end

ksu = ny*(nx+1);
ksv = (ny+1)*nx;
for i = 1:nx
  for j = 1:ny  
    for k = 2:nz+1
      udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i-1;  
      vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i-1 + ndu;
      wdx = (k-1)*ny*nx     + (j-1)*nx     + i-1 + ndu +ndv;
      
      %% wf = -zx*um+-zy*vm + (1+zx^2 +zy^2)*wm
       if k<nz+1
        T(wdx+1,wdx+1)= 1 + 0.5*(zx(i,j,k-1)^2 + zy(i,j,k-1)^2 +...
                                 zx(i,j,k  )^2 + zy(i,j,k  )^2 );

        T(wdx+1,udx+1   -ksu)=  -0.25*zxdy(i,j,k-1)*dxu(i  ,j)/Aw(i,j);
        T(wdx+1,udx+1       )=  -0.25*zxdy(i,j,k  )*dxu(i  ,j)/Aw(i,j);
        T(wdx+1,udx+1+is-ksu)=  -0.25*zxdy(i,j,k-1)*dxu(i+1,j)/Aw(i,j);
        T(wdx+1,udx+1+is    )=  -0.25*zxdy(i,j,k  )*dxu(i+1,j)/Aw(i,j);
        
        T(wdx+1,vdx+1   -ksv)=  -0.25*zydx(i,j,k-1)*dyv(i,j  )/Aw(i,j);
        T(wdx+1,vdx+1       )=  -0.25*zydx(i,j,k  )*dyv(i,j  )/Aw(i,j);
        T(wdx+1,vdx+1+js-ksv)=  -0.25*zydx(i,j,k-1)*dyv(i,j+1)/Aw(i,j);
        T(wdx+1,vdx+1+js    )=  -0.25*zydx(i,j,k  )*dyv(i,j+1)/Aw(i,j);
       else  
        T(wdx+1,wdx+1)= 1 + zx(i,j,k-1)^2 + zy(i,j,k-1)^2;
  
        T(wdx+1,udx+1   -ksu)=  -0.5*zxdy(i,j,k-1)*dxu(i  ,j)/Aw(i,j);
        T(wdx+1,udx+1+is-ksu)=  -0.5*zxdy(i,j,k-1)*dxu(i+1,j)/Aw(i,j);
        
        T(wdx+1,vdx+1   -ksv)=  -0.5*zydx(i,j,k-1)*dyv(i,j  )/Aw(i,j);
        T(wdx+1,vdx+1+js-ksv)=  -0.5*zydx(i,j,k-1)*dyv(i,j+1)/Aw(i,j);
          
      end

    end      
  end
end
toc
else
   disp('Fast makeT')
   tic 

   %% Coefficients for Uf   
   cuwcu = zeros(nx+1,ny,nz);
   cuwwu = zeros(nx+1,ny,nz);
   cuwcc = zeros(nx+1,ny,nz);
   cuwwc = zeros(nx+1,ny,nz);
   
   cuucc = ones(nx+1,ny,nz);
   cuucc(nx+1,:,:) = 0;  %%  It doesn't matter because D does not reach nx+1 
   cuwwu(2:nx+1,:,:) = -0.25.*zxdy.*dzw(:,:,2:nz+1)./Au(2:nx+1,:,:);
   cuwcu(1:nx  ,:,:) = -0.25.*zxdy.*dzw(:,:,2:nz+1)./Au(1:nx  ,:,:);
   cuwwc(2:nx+1,:,:) = -0.25.*zxdy.*dzw(:,:,1:nz  )./Au(2:nx+1,:,:);   
   cuwcc(1:nx  ,:,:) = -0.25.*zxdy.*dzw(:,:,1:nz  )./Au(1:nx  ,:,:); 
   
   %% top bc
   cuwwu(:,:,nz) = 2*cuwwu(:,:,nz);
   cuwcu(:,:,nz) = 2*cuwcu(:,:,nz);
   
   %% bottom bc
    

   
   cuwwc(:,:,1) = 0.0;
   cuwcc(:,:,1) = 0.0;
   
   cuvcc = zeros(nx+1,ny);
   cuvwc = zeros(nx+1,ny);
   cuvcn = zeros(nx+1,ny);
   cuvwn = zeros(nx+1,ny);
   
   if strcmp(rname,'FMbc')
     disp('Yes bottom BC')
     if ny>2&nx>2
       beta = 1.0 - 0.5*zx(:,:,1).^2./alpha(:,:,1) ;
       cuucc(2:nx,:,1) = 0.5*(beta(1:nx-1,:) + beta(2:nx,:));
       
       cuvcc(1:nx  ,:) = -zxzy8dzcwi.*dyv(:,1:ny  )./Au(1:nx  ,:,1);
       cuvwc(2:nx+1,:) = -zxzy8dzcwi.*dyv(:,1:ny  )./Au(2:nx+1,:,1);
       cuvcn(1:nx  ,:) = -zxzy8dzcwi.*dyv(:,2:ny+1)./Au(1:nx  ,:,1);
       cuvwn(2:nx+1,:) = -zxzy8dzcwi.*dyv(:,2:ny+1)./Au(2:nx+1,:,1);
     end
   else
     disp('NO bottom BC')
     cuucc(2:nx,:,1) = 1.0;
   end
   
   %% Coefficients for Vf   
   cvwcu = zeros(nx,ny+1,nz);
   cvwsu = zeros(nx,ny+1,nz);
   cvwcc = zeros(nx,ny+1,nz);
   cvwsc = zeros(nx,ny+1,nz);
   
   cvvcc = ones(nx,ny+1,nz);
   cvvcc(:,ny+1,:) = 0;   %%  It doesn't matter because D does not reach ny+1
   cvwsu(:,2:ny+1,:) = -0.25.*zydx.*dzw(:,:,2:nz+1)./Av(:,2:ny+1,:);
   cvwcu(:,1:ny  ,:) = -0.25.*zydx.*dzw(:,:,2:nz+1)./Av(:,1:ny  ,:);
   cvwsc(:,2:ny+1,:) = -0.25.*zydx.*dzw(:,:,1:nz  )./Av(:,2:ny+1,:);   
   cvwcc(:,1:ny  ,:) = -0.25.*zydx.*dzw(:,:,1:nz  )./Av(:,1:ny  ,:); 
  
   %% top bc
   cvwsu(:,:,nz) = 2*cvwsu(:,:,nz);
   cvwcu(:,:,nz) = 2*cvwcu(:,:,nz);
   
   %% bottom bc
    
   beta = 1.0 - 0.5*zy(:,:,1).^2./alpha(:,:,1) ;
   cvvcc(:,2:ny,1) = 0.5*(beta(:,1:ny-1) + beta(:,2:ny));
   
   cvwsc(:,:,1) = 0.0;
   cvwcc(:,:,1) = 0.0;
   
   cvucc = zeros(nx,ny+1);
   cvuec = zeros(nx,ny+1);
   cvucs = zeros(nx,ny+1);
   cvues = zeros(nx,ny+1);
   
   if strcmp(rname,'FMbc')
     if ny>2&nx>2
       cvucc(:,1:ny  ) = -zxzy8dzcwi.*dxu(1:nx  ,:)./Av(:,1:ny  ,1);
       cvuec(:,1:ny  ) = -zxzy8dzcwi.*dxu(2:nx+1,:)./Av(:,1:ny  ,1);
       cvucs(:,2:ny+1) = -zxzy8dzcwi.*dxu(1:nx  ,:)./Av(:,2:ny+1,1);
       cvues(:,2:ny+1) = -zxzy8dzcwi.*dxu(2:nx+1,:)./Av(:,2:ny+1,1);
     end
   else
     cvvcc(:,2:ny,1) = 1.0;
   end
   
   %% Coefficients for Wf   

   cwwcc = zeros(nx,ny,nz+1);
   cwwcc(:,:,2:nz) = 0.5*(alpha(:,:,1:nz-1)+alpha(:,:,2:nz));
   cwwcc(:,:,nz+1) = alpha(:,:,nz);  
   
   cwucc = zeros(nx,ny,nz+1);
   cwuec = zeros(nx,ny,nz+1);
   cwucd = zeros(nx,ny,nz+1);
   cwued = zeros(nx,ny,nz+1);
   
   cwucc(:,:,1:nz  ) = -0.25*zxdy.*repmat(dxu(1:nx  ,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwuec(:,:,1:nz  ) = -0.25*zxdy.*repmat(dxu(2:nx+1,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwucd(:,:,2:nz+1) = -0.25*zxdy.*repmat(dxu(1:nx  ,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwued(:,:,2:nz+1) = -0.25*zxdy.*repmat(dxu(2:nx+1,:),[1 1 nz])./repmat(Aw,[1 1 nz]);
   
   cwvcc = zeros(nx,ny,nz+1);
   cwvnc = zeros(nx,ny,nz+1);
   cwvcd = zeros(nx,ny,nz+1);
   cwvnd = zeros(nx,ny,nz+1);
   
   cwvcc(:,:,1:nz  ) = -0.25*zydx.*repmat(dyv(:,1:ny  ),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvnc(:,:,1:nz  ) = -0.25*zydx.*repmat(dyv(:,2:ny+1),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvcd(:,:,2:nz+1) = -0.25*zydx.*repmat(dyv(:,1:ny  ),[1 1 nz])./repmat(Aw,[1 1 nz]);
   cwvnd(:,:,2:nz+1) = -0.25*zydx.*repmat(dyv(:,2:ny+1),[1 1 nz])./repmat(Aw,[1 1 nz]); 
   
   cwwcc(:,:,1) = 0;
   cwucc(:,:,1) = 0;cwuec(:,:,1) = 0;cwucd(:,:,1) = 0;cwued(:,:,1) = 0;
   cwvcc(:,:,1) = 0;cwvnc(:,:,1) = 0;cwvcd(:,:,1) = 0;cwvnd(:,:,1) = 0;   
   
   cwucd(:,:,nz+1) = 2*cwucd(:,:,nz+1);
   cwued(:,:,nz+1) = 2*cwued(:,:,nz+1);
   cwvcd(:,:,nz+1) = 2*cwvcd(:,:,nz+1);
   cwvnd(:,:,nz+1) = 2*cwvnd(:,:,nz+1);    
           
   
   
   %% Put coefficients in matrix   
   nnz = 5*ndu + 5*ndv +9*ndw +2*4*nx*ny;
   iT = zeros(nnz,1);
   jT = zeros(nnz,1);
   cT = zeros(nnz,1); 
   v = 1;
   
   is = 1;
   js = nx;
   ks = nx*ny;
       
   %% Coefficients for Uf
   [i,j,k] = meshgrid([1:nx+1],[1:ny],[1:nz]);
   i = reshape(permute(i,[2 1 3]),ndu,1);
   j = reshape(permute(j,[2 1 3]),ndu,1);
   k = reshape(permute(k,[2 1 3]),ndu,1);   
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
   wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
  

   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuucc,udx,udx);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuwcc,udx,wdx   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuwwc,udx,wdx-is);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuwcu,udx,wdx   +ks);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuwwu,udx,wdx-is+ks);
    
   %% extra bottom cross terms
   jsv = nx;
   ndbu = ny*(nx+1);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuvwc,udx(1:ndbu),vdx(1:ndbu)-is);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuvcc,udx(1:ndbu),vdx(1:ndbu)   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuvwn,udx(1:ndbu),vdx(1:ndbu)-is+jsv);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cuvcn,udx(1:ndbu),vdx(1:ndbu)   +jsv);


   %% Coefficients for Vf
   [i,j,k] = meshgrid([1:nx],[1:ny+1],[1:nz]);
   i = reshape(permute(i,[2 1 3]),ndv,1);
   j = reshape(permute(j,[2 1 3]),ndv,1);
   k = reshape(permute(k,[2 1 3]),ndv,1);   
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
   wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
  
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvvcc,vdx,vdx);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvwcc,vdx,wdx   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvwsc,vdx,wdx-js);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvwcu,vdx,wdx   +ks);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvwsu,vdx,wdx-js+ks);
    
   %% extra bottom cross terms
   jsu = nx+1;
   ndbv = (ny+1)*nx;
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvucc,vdx(1:ndbv),udx(1:ndbv)   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvuec,vdx(1:ndbv),udx(1:ndbv)+is);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvucs,vdx(1:ndbv),udx(1:ndbv)   -jsu);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cvues,vdx(1:ndbv),udx(1:ndbv)+is-jsu);  
   
   %% Coefficients for Wf
   [i,j,k] = meshgrid([1:nx],[1:ny],[1:nz+1]);
   i = reshape(permute(i,[2 1 3]),ndw,1);
   j = reshape(permute(j,[2 1 3]),ndw,1);
   k = reshape(permute(k,[2 1 3]),ndw,1);   
   udx = (k-1)*ny*(nx+1) + (j-1)*(nx+1) + i;  
   vdx = (k-1)*(ny+1)*nx + (j-1)*nx     + i + ndu;
   wdx = (k-1)*ny*nx     + (j-1)*nx     + i + ndu +ndv;
 
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwwcc,wdx,wdx   );
   
   ksu = ny*(nx+1);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwucc,wdx,udx   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwuec,wdx,udx+is);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwucd,wdx,udx   -ksu);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwued,wdx,udx+is-ksu);  
   
   ksv = (ny+1)*nx;
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwvcc,wdx,vdx   );
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwvnc,wdx,vdx+js);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwvcd,wdx,vdx   -ksu);
   [iT,jT,cT,v] = addcoef(iT,jT,cT,v,cwvnd,wdx,vdx+js-ksu); 
    
   
   iT(cT==0) = [];
   jT(cT==0) = [];
   cT(cT==0) = [];
   nnz = length(cT);  
   
   T = sparse(iT,jT,cT,ndu+ndv+ndw,ndu+ndv+ndw,nnz);
toc
end
