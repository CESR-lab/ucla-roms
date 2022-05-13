% fill_frc_pipes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Put in the daily Pipe runoff data in the 'new' format. 
%     
%
%  2020, Jeroen Molemaker (UCLA)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%

if 1
 clear all
 close all
 wrk_dir = '/avatar/nmolem/WEC/';
 grdname   = [wrk_dir, 'sample_grd.nc'];
 frcname   = [wrk_dir, 'sample_wwv.nc'];
 cdistname = [wrk_dir, 'sample_cdist.mat'];
end
%
%  Pipe flux data file names:
%
pipes_dir = '/avatar/nmolem/batavia/nmolem/OBSERV/RUNOFF/';
%runoff_dir = './';
pipes_data  = [pipes_dir,'pipes_gom_2010.nc'];

mask = ncread(grdname,'mask_rho');
[nx ny] = size(mask);


 if 0 % create entry for pipe distribution in the grid file if needed
   nccreate(grdname,'pipe_flux','Dimensions',{'xi_rho',nx,'eta_rho',ny},'datatype','single');
   ncwriteatt(grdname,'pipe_flux','long_name','pipe volume flux partition');
   ncwriteatt(grdname,'pipe_flux','units','none');
   rflx = 0*mask;
   ncwrite(grdname,'pipe_flux',rflx);
 end



%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% disp(' ')
% disp(' Read in the grid_data')
 lon  = ncread(grdname,'lon_rho');
 lat  = ncread(grdname,'lat_rho');
 dxi  = ncread(grdname,'pm');
 dyi  = ncread(grdname,'pn');
 mask = ncread(grdname,'mask_rho');
 h    = ncread(grdname,'hraw');
 dx  = 1./dxi;
 dy  = 1./dyi;

 [nx, ny] = size(lon);
 

 if 0 % reset pipe flux if needed
   rflx = 0*mask;
   ncwrite(grdname,'pipe_flux',rflx);
 end

% find the appropriate location of the pipe in the grid.

if 0
 lon_frc = ncread(pipes_data,'lon');
 lat_frc = ncread(pipes_data,'lat');
 tim_frc = ncread(pipes_data,'rno_time');
 flx_frc = ncread(pipes_data,'runoff');
 [nt,nr] = size(flx_frc);
 else
  % create your own pipes
  nt = 360; % once per day
  nr =   1; % 1 is plenty
  lon_frc(1) = -120.670;
  lat_frc(1) =  34.515;
  tim_frc = [0.5:1:359.5]; % time in days
  flx_frc = 2e3*ones(360,1); % m3/s

 end


 d2r = pi/180;
 iloc = 0*lon_frc;
 jloc = 0*lon_frc;
 mindist = 0*lon_frc;
 for i = 1:length(lon_frc)
   dist = gc_dist(lon*d2r,lat*d2r,lon_frc(i)*d2r,lat_frc(i)*d2r);
   mindist(i) = min(dist(:));
   if mindist(i) < mean(dx(:))
     [iloc(i),jloc(i)] = find(dist==mindist(i));
   end
 end



 % toss the pipes outside the grid
 out = mindist>mean(dx(:));
 iloc(out) = [];
 jloc(out) = [];
 mindist(out) = [];
 npip = length(iloc) 
 pip_flx = zeros(nt,npip);
 pip_lon = zeros(npip);
 pip_lat = zeros(npip);
 ipip = 0;
 for i = 1:nr
   if ~out(i)
     ipip = ipip+1;
     pip_flx(:,ipip) = flx_frc(:,i);
     pip_lon(ipip)   = lon_frc(i);
     pip_lat(ipip)   = lat_frc(i);
   end
 end
 pip_tim = tim_frc;

% Edit the 2d field with pipe fluxes


 for ipip = 1:npip

	 ipip
   i0 = max(iloc(ipip)-50,1);
   i1 = min(iloc(ipip)+50,nx);
   j0 = max(jloc(ipip)-50,1);
   j1 = min(jloc(ipip)+50,ny);
   if ipip==3
     i0 = i0-600;
     i1 = i1+800;
     j0 = j0-200;
     j1 = j1+200;j1 = min(j1,ny)
   end

   % edit the mask if needed

   % edit the pipe flux
   pipe_click
 end

 % Check the 2d field with pipe fluxes and partition the flux
 mask = ncread(grdname,'mask_rho');
 rflx = ncread(grdname,'pipe_flux');
 for ipip = 1:npip
   ipip
   err_rflx = rflx>0&mask<1;
   if sum(err_rflx(:))
     rflx(err_rflx) = 0;
     disp('some fluxes were wrong and were removed')
   end

   sum_flx = sum(rflx(:)>=ipip&rflx(:)<ipip+1) 
   if sum_flx>0
     rflx(rflx>=ipip&rflx<ipip+1) = ipip + 1./sum_flx;
   end
 end
 ncwrite(grdname,'pipe_flux',rflx);


 % Done with the pipe partitions, now write the volumes and tracer values


 ntimes = 360;
 nt = 2; % Temperature and Salinity only

 pip_trc = zeros(npip,nt,ntimes);
 for ipip = 1:npip
%  pip_trc(ipip,1,:) = 5 + 13*(cos(tim_frc*2*pi/365 + 0.85*pi)+1); % idealized annual cycle of temp
   pip_trc(ipip,1,:) = 17;
   pip_trc(ipip,2,:) = 1.0; % Salinity
 end

 if 0 % create entry for pipe volume and tracer data in the forcing file if needed
   nccreate(frcname,'pipe_volume','Dimensions',{'npipe',npip,'pipe_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_volume','long_name','River volume flux');
   ncwriteatt(frcname,'pipe_volume','units','m^3/s');

   nccreate(frcname,'pipe_tracer','Dimensions',{'npipe',npip,'ntracers',nt,'pipe_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_tracer','long_name','River tracer data');
   ncwriteatt(frcname,'pipe_tracer','units','variable');

   nccreate(frcname,'pipe_time','Dimensions',{'pipe_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_time','long_name','pipe data time');
   ncwriteatt(frcname,'pipe_time','units','yearday');
 end

 cycle = 1;  %% cyclical
 if (cycle)
   ncwriteatt(frcname,'pipe_time','cycle_length',360);
 end

 ncwrite(frcname,'pipe_volume',pip_flx');
 ncwrite(frcname,'pipe_tracer',pip_trc);
 ncwrite(frcname,'pipe_time',pip_tim');

