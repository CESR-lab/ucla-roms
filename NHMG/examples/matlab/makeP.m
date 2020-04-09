%% construct P directly

 %   ! Define matrix coefficients cA
 %   ! Coefficients are stored in order of diagonals (Fortran code)  
 %   ! cA(1,:,:,:)      -> p(k,j,i)
 %   ! cA(2,:,:,:)      -> p(k-1,j,i)
 %   ! cA(3,:,:,:)      -> p(k+1,j-1,i)
 %   ! cA(4,:,:,:)      -> p(k,j-1,i)
 %   ! cA(5,:,:,:)      -> p(k-1,j-1,i)
 %   ! cA(6,:,:,:)      -> p(k+1,j,i-1)
 %   ! cA(7,:,:,:)      -> p(k,j,i-1)
 %   ! cA(8,:,:,:)      -> p(k-1,j,i-1)
 
    alpha = (1 + zx.^2 + zy.^2); %% defined at rho-points
    zxzy8dzcwi = 0.125*zx(:,:,1).*zy(:,:,1).*dz(:,:,1)./alpha(:,:,1);

    cA = zeros(8,nx,ny,nz);

    for i = 1:nx   
      for j = 1:ny
        for k = 2:nz
             cA(2,i,j,k) =  Aw(i,j)/dzw(i,j,k)*(1+ ...                                 %%  couples with k-1
                  0.5*(zx(i,j,k-1)^2+zy(i,j,k-1)^2+zx(i,j,k  )^2+zy(i,j,k  )^2 ));
        end
      end
    end
    for i = 1:nx
      for j = 2:ny
        for k = 2:nz-1 
             cA(3,i,j,k) =  0.25*( zy(i,j,k+1)*dx(i,j)+ zy(i,j-1,k)*dx(i,j-1) );       %% couples with k+1 j-1
             cA(4,i,j,k) =  Av(i,j,k)/dyv(i,j);                                        %% couples with j-1
             cA(5,i,j,k) = -0.25*( zy(i,j,k-1)*dx(i,j)+ zy(i,j-1,k)*dx(i,j-1) );       %% couples with k-1 j-1
        end
      end
    end
    
    for i = 2:nx
      for j = 1:ny
        for k = 2:nz-1 
          cA(6,i,j,k) =  0.25*( zx(i,j,k+1)*dy(i,j)+ zx(i-1,j,k)*dy(i-1,j) );       %% couples with k+1 i-1
          cA(7,i,j,k) =  Au(i,j,k)/dxu(i,j);                                        %% couples with i-1
          cA(8,i,j,k) = -0.25*( zx(i,j,k-1)*dy(i,j)+ zx(i-1,j,k)*dy(i-1,j) );       %% couples with k-1 i-1
        end
      end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k = 1;  %% bottom level
    for i = 1:nx
      for j = 2:ny
         cA(3,i,j,k) =  0.25*( zy(i,j,k+1)*dx(i,j)+ zy(i,j-1,k)*dx(i,j-1) );       %% couples with k+1 j-1
         
         alp = 1-0.25*(zy(i,j,1)^2/alpha(i,j,1)+zy(i,j-1,1)^2/alpha(i,j-1,1)); 
         cA(4,i,j,k) =  alp*Av(i,j,k)/dyv(i,j) ...                                 %% couples with j-1    
                        - 0.25*zy(i,j-1,k)*dx(i,j-1)+0.25*zy(i,j,k)*dx(i,j);
      end
    end
    for i = 2:nx
      for j = 1:ny
         cA(6,i,j,k) =  0.25*( zx(i,j,k+1)*dy(i,j)+ zx(i-1,j,k)*dy(i-1,j) );       %% couples with k+1 i-1 
        
         alp = 1-0.25*(zx(i,j,1)^2/alpha(i,j,1)+zx(i-1,j,1)^2/alpha(i-1,j,1)); 
         cA(7,i,j,k) =  alp*Au(i,j,k)/dxu(i,j) ...                                 %% couples with i-1    
                       - 0.25*zx(i-1,j,k)*dy(i-1,j)+0.25*zx(i,j,k)*dy(i,j); 
      end
    end 
    
    for i = 2:nx
      for j = 1:ny-1
         cA(5,i,j,k) =   zxzy8dzcwi(i-1,j,k) + zxzy8dzcwi(i,j+1,k);              %% couples with j+1 i-1 
      end
    end    
    for i = 2:nx
      for j = 2:ny
         cA(8,i,j,k) =  -zxzy8dzcwi(i-1,j,k) - zxzy8dzcwi(i,j-1,k);              %% couples with j-1 i-1
      end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k = nz;  %% top level
    for i = 1:nx
      for j = 2:ny
          
         cA(4,i,j,k) =  Av(i,j,k)/dyv(i,j) ...                                     %% couples with j-1    
                        - 0.25*zy(i,j-1,k)*dx(i,j-1)+0.25*zy(i,j,k)*dx(i,j);
         cA(5,i,j,k) = -0.25*( zy(i,j,k-1)*dx(i,j)+ zy(i,j-1,k)*dx(i,j-1) );       %% couples with k-1 j-1
      end
    end    
    for i = 2:nx
      for j = 1:ny         
         cA(7,i,j,k) =  Au(i,j,k)/dxu(i,j) ...                                     %% couples with i-1    
                        - 0.25*zx(i-1,j,k)*dy(i-1,j)+0.25*zx(i,j,k)*dy(i,j);
         cA(8,i,j,k) = -0.25*( zx(i,j,k-1)*dy(i,j)+ zx(i-1,j,k)*dy(i-1,j) );       %% couples with k-1 i-1
      end
    end    
return

%%% Top Boundary:  remember that dzw(nz+1) is half a normal dzw

%           u+: ( p(i+1,j,k  ) - p(i  ,j,k  ))*Au(i+1,j,k)/dxu(i+1,j)
%                              + p(i  ,j,k  )*2*zx4dy(i  ,j,k) ...
%                              + p(i+1,j,k  )*2*zx4dy(i+1,j,k) ...
%              -( p(i  ,j,k  ) - p(i  ,j,k-1))*zx4dy(i  ,j,k) ...
%              -( p(i+1,j,k  ) - p(i+1,j,k-1))*zx4dy(i+1,j,k) ...

%           u-: ( p(i  ,j,k  ) - p(i-1,j,k  ))*alp*Au(i,j,k)/dxu(i,j)
%                              + p(i-1,j,k  )*2*zx4dy(i-1,j,k) ...
%                              + p(i  ,j,k  )*2*zx4dy(i  ,j,k) ...
%              -( p(i-1,j,k+1) - p(i-1,j,k  ))*zx4dy(i-1,j,k) ...
%              -( p(i  ,j,k+1) - p(i  ,j,k  ))*zx4dy(i  ,j,k) ...

%           v+: ( p(i,j+1,k  ) - p(i,j  ,k  ))*alp*Av(i,j+1,k)/dyv(i,j+1)
%                              + p(i,j  ,k  )*2*zy4dx(i,j  ,k) ...
%                              + p(i,j+1,k  )*2*zy4dx(i,j+1,k) ...
%              -( p(i,j  ,k  ) - p(i,j  ,k-1))*zy4dx(i,j  ,k) ...
%              -( p(i,j+1,k  ) - p(i,j+1,k-1))*zy4dx(i,j+1,k) ...
%
%           v-: ( p(i,j  ,k  ) - p(i,j-1,k  ))*Av(i,j,k)/dyv(i,j)
%                              + p(i,j-1,k  )*2*zy4dx(i,j-1,k) ...
%                              + p(i,j  ,k  )*2*zy4dx(i,j  ,k) ...
%              -( p(i,j-1,k  ) - p(i,j-1,k-1))*zy4dx(i,j-1,k) ...
%              -( p(i,j  ,k  ) - p(i,j  ,k-1))*zy4dx(i,j  ,k) ...
%  
%           alp =  1 + 0.5*(zx(i,j,k+1)^2 + zy(i,j,k+1)^2 +...
%                                 zx(i,j,k  )^2 + zy(i,j,k  )^2 );
%           w+: (     - p(i  ,j,k  ))*alp*Aw(i,j,k+1)/dzwi(i,j,k+1) ...
%              -( p(i  ,j,k  ) - p(i-1,j,k  ))*2*zx4dy(i,j,k  ) ...
%              -( p(i+1,j,k  ) - p(i  ,j,k  ))*2*zx4dy(i,j,k  ) ...
%              -( p(i,j  ,k  ) - p(i,j-1,k  ))*2*zy4dx(i,j,k  ) ...
%              -( p(i,j+1,k  ) - p(i,j  ,k  ))*2*zy4dx(i,j,k  ) ...
%
%           alp =  1 + 0.5*(zx(i,j,k-1)^2 + zy(i,j,k-1)^2 +...
%                                 zx(i,j,k  )^2 + zy(i,j,k  )^2 );
%           w-: ( p(i  ,k  ) - p(i  ,k-1))*alp*Aw(i,j,k)/dzwi(i,j,k) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4dy(i,j,k  ) ...
%              -( p(i  ,k-1) - p(i-1,k-1))*zx4dy(i,j,k-1) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4dy(i,j,k  ) ...
%              -( p(i+1,k-1) - p(i  ,k-1))*zx4dy(i,j,k-1) ...

%              -( p(i,j  ,k  ) - p(i,j-1,k  ))*zy4dx(i,j,k  ) ...
%              -( p(i,j  ,k-1) - p(i,j-1,k-1))*zy4dx(i,j,k-1) ...
%              -( p(i,j+1,k  ) - p(i,j  ,k  ))*zy4dx(i,j,k  ) ...
%              -( p(i,j+1,k-1) - p(i,j  ,k-1))*zy4dx(i,j,k-1);
%




%%% Bottom Boundary:

%           alp = (1-0.25*(zx(i,j,1)^2*cwi(i,j,k)+zx(i+1,j,1)^2*cwi(i+1,j,k))
%           u+: ( p(i+1,j,k  ) - p(i  ,j,k  ))*alp*Au(i+1,j,k)/dxu(i+1,j)
%              -( p(i  ,j,k+1) - p(i  ,j,k  ))*zx4dy(i  ,j,k) ...
%              -( p(i+1,j,k+1) - p(i+1,j,k  ))*zx4dy(i+1,j,k) ...
%              -( p(i+1,j+1,k) - p(i+1,j  ,k))*zxzy8dzcwi(i+1,j,k) ...
%              -( p(i+1,j  ,k) - p(i+1,j-1,k))*zxzy8dzcwi(i+1,j,k) ...
%              -( p(i  ,j+1,k) - p(i  ,j  ,k))*zxzy8dzcwi(i  ,j,k) ...
%              -( p(i  ,j  ,k) - p(i  ,j-1,k))*zxzy8dzcwi(i  ,j,k) ...

%           alp = (1-0.25*(zx(i,j,1)^2*cwi(i,j,k)+zx(i-1,j,1)^2*cwi(i-1,j,k))
%           u-: ( p(i  ,j,k  ) - p(i-1,j,k  ))*alp*Au(i,j,k)/dxu(i,j)
%              -( p(i-1,j,k+1) - p(i-1,j,k  ))*zx4dy(i-1,j,k) ...
%              -( p(i  ,j,k+1) - p(i  ,j,k  ))*zx4dy(i  ,j,k) ...
%              -( p(i  ,j+1,k) - p(i  ,j  ,k))*zxzy8dzcwi(i  ,j,k) ...
%              -( p(i  ,j  ,k) - p(i  ,j-1,k))*zxzy8dzcwi(i  ,j,k) ...
%              -( p(i-1,j+1,k) - p(i-1,j  ,k))*zxzy8dzcwi(i-1,j,k) ...
%              -( p(i-1,j  ,k) - p(i-1,j-1,k))*zxzy8dzcwi(i-1,j,k) ...

%           alp = (1-0.25*(zy(i,j,1)^2*cwi(i,j,k)+zy(i,j+1,1)^2*cwi(i,j+1,k))
%           v+: ( p(i,j+1,k  ) - p(i,j  ,k  ))*alp*Av(i,j+1,k)/dyv(i,j+1)
%              -( p(i,j  ,k+1) - p(i,j  ,k  ))*zy4dx(i,j  ,k) ...
%              -( p(i,j+1,k+1) - p(i,j+1,k  ))*zy4dx(i,j+1,k) ...
%              -( p(i+1,j+1,k) - p(i  ,j+1,k))*zxzy8dzcwi(i,j+1,k) ...
%              -( p(i  ,j+1,k) - p(i-1,j+1,k))*zxzy8dzcwi(i,j+1,k) ...
%              -( p(i+1,j  ,k) - p(i  ,j  ,k))*zxzy8dzcwi(i,j  ,k) ...
%              -( p(i  ,j  ,k) - p(i-1,j  ,k))*zxzy8dzcwi(i,j  ,k) ...
%
%           alp = (1-0.25*(zy(i,j,1)^2*cwi(i,j,k)+zy(i,j-1,1)^2*cwi(i,j-1,k))
%           v-: ( p(i,j  ,k  ) - p(i,j-1,k  ))*alp*Av(i,j,k)/dyv(i,j)
%              -( p(i,j-1,k+1) - p(i,j-1,k  ))*zy4dx(i,j-1,k) ...
%              -( p(i,j  ,k+1) - p(i,j  ,k  ))*zy4dx(i,j  ,k) ...
%              -( p(i+1,j  ,k) - p(i  ,j  ,k))*zxzy8dzcwi(i,j  ,k) ...
%              -( p(i  ,j  ,k) - p(i-1,j  ,k))*zxzy8dzcwi(i,j  ,k) ...
%              -( p(i+1,j-1,k) - p(i  ,j-1,k))*zxzy8dzcwi(i,j-1,k) ...
%              -( p(i  ,j-1,k) - p(i-1,j-1,k))*zxzy8dzcwi(i,j-1,k) ...
%
%  
%           alp =  1 + 0.5*(zx(i,j,k+1)^2 + zy(i,j,k+1)^2 +...
%                                 zx(i,j,k  )^2 + zy(i,j,k  )^2 );
%           w+: ( p(i  ,k+1) - p(i  ,k  ))*alp*Aw(i,j,k+1)/dzwi(i,j,k+1) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4dy(i,j,k  ) ...
%              -( p(i  ,k+1) - p(i-1,k+1))*zx4dy(i,j,k+1) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4dy(i,j,k  ) ...
%              -( p(i+1,k+1) - p(i  ,k+1))*zx4dy(i,j,k+1) ...

%              -( p(i,j  ,k  ) - p(i,j-1,k  ))*zy4dx(i,j,k  ) ...
%              -( p(i,j  ,k+1) - p(i,j-1,k+1))*zy4dx(i,j,k+1) ...
%              -( p(i,j+1,k  ) - p(i,j  ,k  ))*zy4dx(i,j,k  ) ...
%              -( p(i,j+1,k+1) - p(i,j  ,k+1))*zy4dx(i,j,k+1);
%
%           w-: 0 






%                           % u1 =  u2 - w2*z_x
%                           % w1 = -zx*u2+(1+zx^2)*w2
%
%    w  u     D*MT*G          u+  - u-  + w+ - w-
%
%
%           u+: ( p(i+1,k  ) - p(i  ,k  ))*Au(i+1,j,k)/dxu(i+1,j)  ...
%              -( p(i  ,k  ) - p(i  ,k-1))*zx4dy(i  ,j,k) ...
%              -( p(i  ,k+1) - p(i  ,k  ))*zx4dy(i  ,j,k) ...
%              -( p(i+1,k  ) - p(i+1,k-1))*zx4dy(i+1,j,k) ... 
%              -( p(i+1,k+1) - p(i+1,k  ))*zx4dy(i+1,j,k);

%
%           u-:-( p(i  ,k  ) - p(i-1,k  ))*Au(i  ,j,k)*dxui(i  ) ...
%              +( p(i-1,k  ) - p(i-1,k-1))*zx4dy(i-1,j,k) ...
%              +( p(i-1,k+1) - p(i-1,k)  )*zx4dy(i-1,j,k) ...
%              +( p(i  ,k  ) - p(i  ,k-1))*zx4dy(i  ,j,k) ...
%              +( p(i  ,k+1) - p(i  ,k)  )*zx4dy(i  ,j,k);
%
%           w+: ( p(i  ,k+1) - p(i  ,k  ))*dx(i)*dzwi(i,k+1)*(1+zxw^2) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              -( p(i  ,k+1) - p(i-1,k+1))*zx4(i,k+1) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  ) ...
%              -( p(i+1,k+1) - p(i  ,k+1))*zx4(i,k+1);
%    
%           wi:-( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i  ,k-1) - p(i-1,k-1))*zx4(i,k-1) ...
%              +( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              +( p(i+1,k-1) - p(i  ,k-1))*zx4(i,k-1) ...
%              +( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  );


ccc = zeros(nx,nz);
for k =1:nz
   ccc(:,k) = - dzu(2:nx+1,k).*dxiu(2:nx+1)-dzu(1:nx,k).*dxiu(1:nx) ...
                - dx.*dziw(:,k+1).*(1+zxw(:,k+1).^2) ...
                - dx.*dziw(:,k  ).*(1+zxw(:,k).^2);
end

ccb = zeros(nx,nz+1);       
for k =1:nz+1
   ccb(:,k) = dx.*dziw(:,k).*(1+zxw(:,k).^2); 
end
cct = ccb(:,2:nz+1);
ccb = ccb(:,1:nz);


cwc = zeros(nx+1,nz);  % for p(i-1,j,k)
for k = 1:nz
   cwc(:,k) =  dzu(:,k).*dxiu;
end
cec = cwc(2:nx+1,:);
cwc = cwc(1:nx,:);

        
cwt = zeros(nx,nz);
cwt(2:nx,:) = 0.25*zx(1:nx-1,:);
cwt(:,1:nz-1) = cwt(:,1:nz-1) + 0.25*zx(:,2:nz);
cwt(:,nz) = cwt(:,nz) + 0.25*zx_sur;

ceb = zeros(nx,nz);
ceb(1:nx-1,:) = 0.25*zx(2:nx,:);
ceb(:,2:nz) = ceb(:,2:nz) + 0.25*zx(:,1:nz-1);
ceb(:,1) = ceb(:,1) + 0.25*zx_bot;

cwb = zeros(nx,nz);
cwb(2:nx,:) = -0.25*zx(1:nx-1,:);
cwb(:,2:nz) = cwb(:,2:nz) - 0.25*zx(:,1:nz-1);
cwb(:,1) = cwb(:,1) - 0.25*zx_bot;

cet = zeros(nx,nz);
cet(1:nx-1,:) = -0.25*zx(2:nx,:);
cet(:,1:nz-1) = cet(:,1:nz-1) - 0.25*zx(:,2:nz);
cet(:,nz) = cet(:,nz) - 0.25*zx_sur;

          
%% bc's at the left:  no flux and no px, slightly subtle!
ccc(1,:) = ccc(1,:) + cwc(1,:);
cct(1,1:nz-1) = cct(1,1:nz-1) + zx4(1,2:nz) - zx4(1,1:nz-1  );
ccb(1,2:nz  ) = ccb(1,2:nz  ) + zx4(1,2:nz) - zx4(1,1:nz-1);

%% bc's at the right
ccc(nx,:) = ccc(nx,:) + cec(nx,:);
cct(nx,1:nz-1) = cct(nx,1:nz-1) - zx4(nx,2:nz) + zx4(nx,1:nz-1);
ccb(nx,2:nz  ) = ccb(nx,2:nz  ) - zx4(nx,2:nz) + zx4(nx,1:nz-1);

%% bc's at the top
ccc(:,nz) = ccc(:,nz) - cct(:,nz);
cwc(:,nz) = cwc(:,nz) - cwt(:,nz);
cec(:,nz) = cec(:,nz) - cet(:,nz);
 
nd = nx*(nz+1);
is = 1;
ks = nx;
A =sparse(nd,nd,9*nd);         
for i = 1:nx
   for k = 1:nz-1
       idx = (k-1)*nx + i-1 +nx;
       if i>1
         A(idx+1,idx+1-is)   = cwc(i,k);
         A(idx+1,idx+1-is+ks)= cwt(i,k);
         A(idx+1,idx+1-is-ks)= cwb(i,k);
       end
       if i<nx
         A(idx+1,idx+1+is)   = cec(i,k);
         A(idx+1,idx+1+is+ks)= cet(i,k);
         A(idx+1,idx+1+is-ks)= ceb(i,k);
       end
       A(idx+1,idx+1-ks)= ccb(i,k);
       A(idx+1,idx+1+ks)= cct(i,k);
       A(idx+1,idx+1   )= ccc(i,k); 
   end
   k = nz;
   idx = (k-1)*nx + i-1 + nx;
   if i>1
      A(idx+1,idx+1-is)   = cwc(i,k);
      A(idx+1,idx+1-is-ks)= cwb(i,k);
   end
   if i<nx
      A(idx+1,idx+1+is)   = cec(i,k);
      A(idx+1,idx+1+is-ks)= ceb(i,k);
   end
   A(idx+1,idx+1-ks)= ccb(i,k);
   A(idx+1,idx+1   )= ccc(i,k);    
end

%% Bottom boundary condition
for i = 1:nx
   idx =  i-1;
   if i>1
     A(idx+1,idx+1-is+ks) =  0.25*zx(i,1);
   end
   if i<nx
      A(idx+1,idx+1+is+ks) = -0.25*zx(i,1);
   end
   A(idx+1,idx+1      ) = -0.5*dx(i)*dziw(i,1)*(1+zxw(i,1)*zxw(i,1));
   
   A(idx+1,idx+1   +ks) = +0.5*dx(i)*dziw(i,1)*(1+zxw(i,1)*zxw(i,1));
   
end

%% Get a coefficient out of the (P4) matrix for comparison
cf = zeros(nx,nz);
is = 1;
ks = nx;
for i = 2:nx-1
   for k = 2:nz-1
      idx = (k-1)*nx + i-1 +nx;
      cf(i,k) = P(idx+1,idx+1+ks);
%      cf(i,k) = P4(idx+1,idx+1+is);   %% couples with p(i,k+1)
   end
end

u = zeros(nx+1,nz);
w = cos(xw*pi/sizex); w(:,1) = 0;
Xr = uw2x(u,w);  %% physical velocity
%Xm = Trm*Xr;    %% physical to momentum (not neccesary here because u==0)
div = D*T*Xr;
Fbot = zeros(nx,1);%
%rhs = [Fbot' div']';
rhs = div;
p = A\rhs;
p2 = reshape(p,[nx nz+1]);

xp = zeros(nx,nz+1);
zp = zeros(nx,nz+1);
xp(:,2:nz+1) = xr;
zp(:,2:nz+1) = zr;
xp(:,1) = xp(:,2);
zp(:,1) = zr_bot;

