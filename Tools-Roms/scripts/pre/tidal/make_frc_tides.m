%----- Begin user input --------------

   tname = 'tpxo9.v2a.nc';
   gname = 'roms_grd.nc';
   fname = 'roms_frc_tides.nc';

   nc = 10; % number of tidal constituents

   dn0 = datenum(1992,1,1);

 % Reference time of simulation
 year = 2000;
 month=    1;
 day  =    1;

%----- End user input --------------

 dn = datenum(year,month,day);

date_mjd=mjd(year,month,day);
[pf,pu,t0,phase_mkB]=egbert_correc(date_mjd,0,0,0);
pu = pu*pi/180;
phase_mkB = phase_mkB*pi/180;
aa = phase_mkB;

t = (dn-dn0)*3600*24;


lon = ncread(gname,'lon_rho');
lat = ncread(gname,'lat_rho');
  h = ncread(gname,'h');
mask= ncread(gname,'mask_rho');
ang = ncread(gname,'angle');
[nx,ny] = size(lon);

 % create forcing file
 om = ncread(tname,'omega',[1],[nc]);
 con= ncreadatt(tname,'/','Constituents');
 if ~exist(fname)
   create_frc_tides(gname,fname,nc,om,con)
   ncwriteatt(fname,'/','Reference_time',datestr(dn));
 end

%  find limits in the tpxo data file
 lon(lon<0) = lon(lon<0)+360;  % tpxo is between 0 and 360;

 lon0 = min(lon(:));
 lon1 = max(lon(:));
 lat0 = min(lat(:));
 lat1 = max(lat(:));

 tlonr = ncread(tname,'lon_r',[1 1],[inf 1]);
 tlatr = ncread(tname,'lat_r',[1 1],[1 inf]);

 i0 = find(tlonr<lon0,1,'last');
 i1 = find(tlonr>lon1,1,'first');
 j0 = find(tlatr<lat0,1,'last');
 j1 = find(tlatr>lat1,1,'first');
 tnx = i1-i0+1;
 tny = j1-j0+1;

 tlonr = ncread(tname,'lon_r',[i0 j0],[tnx 1]);
 tlonu = ncread(tname,'lon_u',[i0 j0],[tnx 1]);
 tlonv = ncread(tname,'lon_v',[i0 j0],[tnx 1]);
 tlatr = ncread(tname,'lat_r',[i0 j0],[1 tny]);
 tlatu = ncread(tname,'lat_u',[i0 j0],[1 tny]);
 tlatv = ncread(tname,'lat_v',[i0 j0],[1 tny]);

 tlon2d = ncread(tname,'lon_r',[i0 j0],[tnx tny]);
 tlat2d = ncread(tname,'lat_r',[i0 j0],[tnx tny]);

 % Read elevations,sal, and transports from tpxo file
 thr = ncread(tname,'h_Re',[i0 j0 1],[tnx tny nc]);
 thi = ncread(tname,'h_Im',[i0 j0 1],[tnx tny nc]);
 thc = complex(thr,thi);
 tur = ncread(tname,'u_Re',[i0 j0 1],[tnx tny nc]);
 tui = ncread(tname,'u_Im',[i0 j0 1],[tnx tny nc]);
 tuc = complex(tur,tui);
 tvr = ncread(tname,'v_Re',[i0 j0 1],[tnx tny nc]);
 tvi = ncread(tname,'v_Im',[i0 j0 1],[tnx tny nc]);
 tvc = complex(tvr,tvi);
 tsr = ncread(tname,'sal_Re',[i0 j0 1],[tnx tny nc]);
 tsi = ncread(tname,'sal_Im',[i0 j0 1],[tnx tny nc]);
 tsc = 1.0*complex(tsr,tsi);   %% The Alan factor is 1.5
 clear 'thr' 'thi' 'tur' 'tui' 'tvr' 'tvi' 'tsr' 'tsi'

 % Get equilibrium tides and correct for SAL
 tpc = equi_tide(tlon2d,tlat2d,nc);
 tpc = tpc - tsc;


 % Multiply with nodal corrections and phase shifts to the reference time
 cI = complex(0,1);
 for ic = 1:nc
  thc(:,:,ic) = pf(ic)*thc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tuc(:,:,ic) = pf(ic)*tuc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tvc(:,:,ic) = pf(ic)*tvc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tpc(:,:,ic) = pf(ic)*tpc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
 end

 % Process variables and write to forcing file
 hc = zeros(nx,ny,nc);
 for ic = 1:nc
   hc(:,:,ic) = interp2(tlonr,tlatr,thc(:,:,ic).',lon,lat);  % non-conjugate transpose: .'
 end
 ncwrite(fname,'ssh_Re',real(hc));
 ncwrite(fname,'ssh_Im',imag(hc));

 % Process tidal barotropic velocities
 uc = zeros(nx,ny,nc);
 vc = zeros(nx,ny,nc);
 for ic = 1:nc
   uc(:,:,ic) = interp2(tlonu,tlatu,tuc(:,:,ic).',lon,lat);
   vc(:,:,ic) = interp2(tlonv,tlatv,tvc(:,:,ic).',lon,lat);
 end
 clear 'tuc' 'tvc'

 % Rotate to grid orientation and convert to barotropic velocity
 cosa = cos(ang);
 sina = sin(ang);
 u = zeros(nx,ny,nc);
 v = zeros(nx,ny,nc);
 for ic = 1:nc
   u(:,:,ic) = (uc(:,:,ic).*cosa + vc(:,:,ic).*sina)./h;
   v(:,:,ic) = (vc(:,:,ic).*cosa - uc(:,:,ic).*sina)./h;
 end

 % Average to u and v points and write to file
 u = 0.5*(u(2:end,:,:)+u(1:end-1,:,:));
 v = 0.5*(v(:,2:end,:)+v(:,1:end-1,:));
 ncwrite(fname,'u_Re',real(u));
 ncwrite(fname,'u_Im',imag(u));
 ncwrite(fname,'v_Re',real(v));
 ncwrite(fname,'v_Im',imag(v));

 pc = zeros(nx,ny,nc);
 for ic = 1:nc
   pc(:,:,ic) = interp2(tlonr,tlatr,tpc(:,:,ic).',lon,lat);
 end
 ncwrite(fname,'pot_Re',real(pc));
 ncwrite(fname,'pot_Im',imag(pc));
return

