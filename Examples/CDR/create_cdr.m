fname = 'cdr_release.nc';
ncdr = 1;
ntracers = 3;

ncid = netcdf.create(fname,'clobber');
dimid = netcdf.defDim(ncid,'xi_rho',1);
dimid = netcdf.defDim(ncid,'eta_rho',1);
netcdf.close(ncid)

% shape of the CDR releases
nccreate(fname,'cdr_lon','Dimensions',{'ncdr',ncdr} );
ncwriteatt(fname,'cdr_lon','long_name','longitude of CDR release');
ncwriteatt(fname,'cdr_lon','units','Degrees East');

nccreate(fname,'cdr_lat','Dimensions',{'ncdr',ncdr} );
ncwriteatt(fname,'cdr_lat','long_name','latitude of CDR release');
ncwriteatt(fname,'cdr_lat','units','Degrees North');

nccreate(fname,'cdr_hsc','Dimensions',{'ncdr',ncdr} );
ncwriteatt(fname,'cdr_hsc','long_name','horizontal scale of CDR release');
ncwriteatt(fname,'cdr_hsc','units','m');

nccreate(fname,'cdr_vsc','Dimensions',{'ncdr',ncdr} );
ncwriteatt(fname,'cdr_vsc','long_name','vertical scale of CDR release');
ncwriteatt(fname,'cdr_vsc','units','m');

nccreate(fname,'cdr_dep','Dimensions',{'ncdr',ncdr} );
ncwriteatt(fname,'cdr_dep','long_name','mean depth of CDR release');
ncwriteatt(fname,'cdr_dep','units','m');

ncwrite(fname,'cdr_lon',-120.718)
ncwrite(fname,'cdr_lat',34.560)
ncwrite(fname,'cdr_hsc',5e3)
ncwrite(fname,'cdr_vsc',10)
ncwrite(fname,'cdr_dep',20)

% Time series of CDR releases
nccreate(fname,'cdr_time','Dimensions',{'time',0} );
ncwriteatt(fname,'cdr_time','long_name','time since 2000-1-1');
ncwriteatt(fname,'cdr_time','units','yearday');
ncwriteatt(fname,'cdr_time','cycle_length',360);

cdr_time = [15:30:345]
ncwrite(fname,'cdr_time',cdr_time);

nccreate(fname,'cdr_volume','Dimensions',{'ncdr',ncdr,'time',0} );
ncwriteatt(fname,'cdr_volume','long_name','Volume flux of CDR release');
ncwriteatt(fname,'cdr_volume','units','m3/s');

cdr_vol = zeros(1,12);
cdr_vol = 0*cdr_vol + 2000;
ncwrite(fname,'cdr_volume',cdr_vol);

nccreate(fname,'cdr_tracer','Dimensions',{'ncdr',ncdr,'ntracers',ntracers,'time',0} );
ncwriteatt(fname,'cdr_tracer','long_name','Tracer concentrations of CDR release');
ncwriteatt(fname,'cdr_tracer','units','C/m3');

cdr_trc = zeros(1,ntracers,12);
cdr_trc(1,1,:) = 20; % temperature
cdr_trc(1,2,:) = 30; % salinity
cdr_trc(1,3,:) =  1; % passive tracer
%ncwrite(fname,'cdr_tracer',cdr_trc);
ncwrite(fname,'cdr_tracer',cdr_trc);

nccreate(fname,'cdr_trcflx','Dimensions',{'ncdr',ncdr,'ntracers',ntracers,'time',0} );
ncwriteatt(fname,'cdr_trcflx','long_name','Tracer fluxes of CDR release');
ncwriteatt(fname,'cdr_trcflx','units','C/s');

% Don't force temp and salt without volume, or you get unintended results
for j = 1:12
  cdr_trc(1,1,j) = 0;
  cdr_trc(1,2,j) = 0;
  cdr_trc(1,3,j) = cdr_vol(1,j)*cdr_trc(1,3,j);
end

ncwrite(fname,'cdr_trcflx',cdr_trc);


