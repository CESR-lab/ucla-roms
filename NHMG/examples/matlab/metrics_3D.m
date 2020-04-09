
if 1
 nx = 50*fac;
 ny = 1;
 nz = 20*fac;
 nx = 25*fac;
 ny = 1*fac;
 nz = 10*fac;
else
  rname = 'FMnbc';
  nx = 10;
  ny = 1; 
  nz = 5;
end

ndu = (nx+1)*ny*nz;
ndv = nx*(ny+1)*nz;
ndw = nx*ny*(nz+1);

%%%input
if 1
  disp('random slopes and grid spacing')
  
  zx = rand(nx,ny,nz);
  zy = rand(nx,ny,nz);

  dx = rand(nx,ny);
  dy = rand(nx,ny);
  dz = rand(nx,ny,nz);
  
  zx(1,:,:)  = 0;
  zx(nx,:,:) = 0;
  zy(:,1,:)  = 0;
  zy(:,ny,:) = 0;
  
else
  disp('2D seamount test case')
  sizex = 1e4;
  sizey = 1e2;
  sizez = 4e3; 
  
  thetas = 1;
  thetab = 0.5;
  
  x1 = sizex*[0.5:1:nx-0.5]'/nx;
  x1u= sizex*[0.0:1:nx-0.0]'/nx;
  y1 = sizey*[0.5:1:ny-0.5]'/ny;
  y1v= sizey*[0.0:1:ny-0.0]'/ny;
  [x,y] = meshgrid(x1,y1);x = x';y=y';
  [xu,~] = meshgrid(x1u,y1);xu = xu';
  [~,yv] = meshgrid(x1,y1v);yv = yv';
  if 1
    x = sizex*(x/sizex+0.1*sin(2*x*pi/sizex) );
    xu= sizex*(xu/sizex+0.1*sin(2*xu*pi/sizex) );
  end
  dx = xu(2:end,:)-xu(1:end-1,:);
  dy = yv(:,2:end)-yv(:,1:end-1); 
  
  xu = repmat(xu,[1 1 nz]);
  xw = repmat(x,[1 1 nz+1]);
  

  
  alpha = 25./(sizex^2);
  h = (sizez - 0.5*sizez*exp(-alpha*(x-0.5*sizex).^2));
  
  
  zr = zlevs(h' ,0*h' ,thetas,thetab,200,nz,'r');
  zw = zlevs(h' ,0*h' ,thetas,thetab,200,nz,'w');
  zr = permute(zr,[3 2 1]);
  zw = permute(zw,[3 2 1]);

  zu = 0*xu;  
  zu(2:nx,:,:) = 0.5*(zr(2:nx,:,:) + zr(1:nx-1,:,:));
  zu(1,:,:) = zu(2,:,:);
  zu(nx+1,:,:) = zu(nx,:,:);
  
  zx = 0*zr;
  zy = 0*zr;
  for k=1:nz
    zx(2:nx-1,:,k) = 0.5*(zr(3:nx,:,k) - zr(1:nx-2,:,k))./dx(2:nx-1,:);
  end
   
  dz = zw(:,:,2:nz+1)-zw(:,:,1:nz);
    
end


%% derived variables

dxu = zeros(nx+1,ny);
dyv = zeros(nx,ny+1);
dzw = zeros(nx,ny,nz+1);

dxu(2:nx,:) = 0.5*(dx(1:nx-1,:) + dx(2:nx,:));
dyv(:,2:ny) = 0.5*(dy(:,1:ny-1) + dy(:,2:ny));
dxu(1,:) = 0.5*dx(1,:);dxu(nx+1,:) = 0.5*dx(nx,:);
dyv(:,1) = 0.5*dy(:,1);dyv(:,ny+1) = 0.5*dy(:,ny);
dzw(:,:,2:nz) =  0.5*(dz(:,:,1:nz-1) + dz(:,:,2:nz));
dzw(:,:,1) = 0.5*dz(:,:,1);dzw(:,:,nz+1) = 0.5*dz(:,:,nz);

zxdy = zeros(nx,ny,nz);
zydx = zeros(nx,ny,nz);
for k = 1:nz
  zxdy(:,:,k) = zx(:,:,k).*dy; 
  zydx(:,:,k) = zy(:,:,k).*dx;
end

Au = zeros(nx+1,ny,nz);
Av = zeros(nx,ny+1,nz);
Aw = zeros(nx,ny);

for k = 1:nz
  Au(2:nx,:,k) = 0.5*( dy(1:nx-1,:).*dz(1:nx-1,:,k)+dy(2:nx,:).*dz(2:nx,:,k) );
  Au(1   ,:,k) = dy(1 ,:).*dz(1 ,:,k);
  Au(nx+1,:,k) = dy(nx,:).*dz(nx,:,k);
  
  Av(:,2:ny,k) = 0.5*( dx(:,1:ny-1).*dz(:,1:ny-1,k)+dx(:,2:ny).*dz(:,2:ny,k) );
  Av(:,1   ,k) = dy(:,1 ).*dz(:,1 ,k);
  Av(:,ny+1,k) = dy(:,ny).*dz(:,ny,k);
  
%  Aw(:,:,k) =  dx.*dy;
end
Aw =  dx.*dy;
