%% Test 3D matrices

metrics_3D % Pre-compute surface areas,lengths and slopes.

% The matrix that transforms momenta into fluxes
T = makeT(zx,zy,zxdy,zydx,Au,Av,Aw,dxu,dyv,dz,dzw,rname);

% The diagonal matrix that has volume elements at u,v, and w points
M = makeM(Au,Av,Aw,dxu,dyv,dzw);

% The inverse of M.  
ndt = ndu+ndv+ndw;
Mi = spdiags(1./spdiags(M),0,ndt,ndt);

% the volume integrated divergence operator
D = makeD(Au,Av,Aw,rname);

if 0  %% Check symmetry of MT (set coefficients to random in metrics_3D)
  MT= M*T;
  MTT = MT';
  Dif1 = MT-MTT;
  Dif1(abs(Dif1)<1e-14) = 0;
end


%% Check symmetry of P
Pt = D*Mi*T'*D';
if 0   %% Check symmetry of P
  Dif2 = Pt-Pt';
  Dif2(abs(Dif2)<1e-14) = 0;
end

%% Carthesian version
Tc = makeTc(zxdy,zydx,dxu,dyv,Aw);
Dc = D*Tc;
Pc = Dc*Mi*Dc';

disp('operators done')

%% Compare constructed P with direct definition of P
makeP

ks = ny*nx;
js = nx;
is = 1;
cAt = zeros(8,nx,ny,nz);
for i = 1:nx
  for j = 1:ny 
    for k = 2:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(1,i,j,k) = Pt(idx,idx   );
      cAt(2,i,j,k) = Pt(idx,idx-ks);
    end
  end
end

for i = 1:nx
  for j = 2:ny 
    for k = 1:nz-1 
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(3,i,j,k) = Pt(idx,idx+ks-js);
    end
  end
end

for i = 1:nx
  for j = 2:ny 
    for k = 1:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(4,i,j,k) = Pt(idx,idx   -js);
    end
  end
end

for i = 1:nx  %% k-1,j-1
  for j = 2:ny 
    for k = 2:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(5,i,j,k) = Pt(idx,idx-ks-js);
    end
  end
end
for i = 2:nx  %% j+1,i-1
  for j = 1:ny-1 
     k=1;
     idx = (k-1)*ny*nx + (j-1)*nx + i;
     cAt(5,i,j,k) = Pt(idx,idx+js-is);
  end
end


for i = 2:nx
  for j = 1:ny 
    for k = 1:nz-1
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(6,i,j,k) = Pt(idx,idx+ks-is);
    end
  end
end
  
for i = 2:nx
  for j = 1:ny 
    for k = 1:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(7,i,j,k) = Pt(idx,idx   -is);
    end
  end
end

for i = 2:nx  %% k-1,i-1
  for j = 1:ny 
    for k = 2:nz
      idx = (k-1)*ny*nx + (j-1)*nx + i;
      cAt(8,i,j,k) = Pt(idx,idx-ks-is);
    end
  end
end
for i = 2:nx  %% j-1,i-1
  for j = 2:ny 
     k=1;
     idx = (k-1)*ny*nx + (j-1)*nx + i;
     cAt(8,i,j,k) = Pt(idx,idx-js-is);
  end
end


cD = cA+cAt;

E2 = max(max(max(abs(cD(2,:,:,:)))))  %% verify Coeff for i,j,k-1
E3 = max(max(max(abs(cD(3,:,:,:)))))  %% verify Coeff for i,j-1,k+1
E4 = max(max(max(abs(cD(4,:,:,:)))))  %% verify Coeff for i,j-1,k
E5 = max(max(max(abs(cD(5,:,:,:)))))  %% verify Coeff for i,j-1,k-1 (plus j-1,i+1)
E6 = max(max(max(abs(cD(6,:,:,:)))))  %% verify Coeff for i-1,j,k+1 
E7 = max(max(max(abs(cD(7,:,:,:)))))  %% verify Coeff for i-1,j,k
E8 = max(max(max(abs(cD(8,:,:,:)))))  %% verify Coeff for i-1,j,k-1 (plus j-1,i-1)


