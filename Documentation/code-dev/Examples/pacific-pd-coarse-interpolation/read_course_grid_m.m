clc ; clear;

% NOTE: roms indices (i,j) = (x,y), but matlab uses matrix indices
%       where (i=row,j=column) and hence (y,x), so need to tranpose!

%% LAT_RHO

% % Read in netcdf coarse grid data.
% file_full_c='/home/devin/code-roms-versions/pdamien/pacsma50km/pacmed_50km_grd.nc';
% nc_full_c_h=ncread(file_full_c,'lat_rho');
% % file_c='/home/devin/code-roms-versions/pdamien/pacsma50km/pacmed_50km_grd.nc';
% % nc_c_h=ncread(file_c,'lat_rho');
%
% % Read in netcdf refined grid data.
% file_full_r='/home/devin/code-roms-versions/pdamien/pacmed25km/pacmed_grd.nc';
% nc_full_r_h=ncread(file_full_r,'lat_rho');
% % file_r='/home/devin/code-roms-versions/pdamien/pacmed25km/pacmed_grd.nc';
% % nc_r_h=ncread(file_r,'lat_rho');
%
% % read in online ROMS interpolation
% file_full_roms='/home/devin/code-roms-versions/pdamien/interpolations/bulk_res_surf_flux.0000.nc-lat-rho';
% nc_full_roms_int=ncread(file_full_roms,'lat_rho');

%% UWND 50km -> 25km

% startloc = [1 1 1];
% count = [ inf inf 1 ];
% % Read in netcdf coarse grid data.
% file_full_c='/home/devin/code-roms-versions/pdamien/pacsma50km/pacmed_50kmriverdis_DPD_uwnd_rename.nc';
% nc_full_c_h=ncread(file_full_c,'uwnd-int',startloc,count);
%
% % Read in netcdf refined grid data.
% file_full_r='/home/devin/code-roms-versions/pdamien/pacmed25km/pacmed_0p25riverdis_corrected_Y2000M01_frc_newQair_uwnd.nc';
% nc_full_r_h=ncread(file_full_r,'uwnd',startloc,count);
%
% % read in online ROMS interpolation
% file_full_roms='/home/devin/code-roms-versions/pdamien/interpolations/bulk_res_surf_flux.0000.nc-uwnd-int-both';
% nc_full_roms_int=ncread(file_full_roms,'uwnd-int',startloc,count);

%% UWND 25km -> 12km

startloc = [1 1 1];
count = [ inf inf 1 ];
% Read in netcdf coarse grid data.
file_full_c='/home/devin/code-roms-versions/pdamien/pacmed25km/pacmed_0p25riverdis_corrected_Y2000M01_frc_newQair_uwnd.nc';
nc_full_c_h=ncread(file_full_c,'uwnd',startloc,count);

% Read in netcdf refined grid data.
file_full_r='/home/devin/code-roms-versions/pdamien/pacbig12km/pacmed_12km_uwnd1.nc';
nc_full_r_h=ncread(file_full_r,'uwnd',startloc,count);

% read in online ROMS interpolation
% file_full_roms='/home/devin/code-roms-versions/pdamien/interpolations/bulk_res_surf_flux.0000.nc-uwnd-int-both';
% nc_full_roms_int=ncread(file_full_roms,'uwnd-int',startloc,count);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Declare variables
% % v_c = zeros( size_c, size_c);                   % coarse variable
% % v_r = zeros( (size_c-1)*2, (size_c-1)*2, 1 );   % Size of refined var_r
% Lmc = 58;
% Mmc = 40;
%
%% Ranges of course grid
% icmin = 0;
% icmax = Lmc-1;                                   % -1 for rho numbering from 0
% jcmin = 0;
% jcmax = Mmc-1;                               % -1 for rho numbering from 0
%
%% create buff variable to simulate read in data from netcdf file
% buff = reshape(nc_c_h,1,[]); % Note no transpose.
%
%% Interpolate coarse grid
% for ic = icmin+1 : icmax-1+1                    % loop to second last course points = ic. % +1 for matlab indices not rho indices
%   for jc = jcmin+1 : jcmax-1+1
%     ir = 2 * ic -1;                             % ir = refined grid point index for equiv. coarse even node, immediate top-right of ic
%     jr = 2 * jc -1;
%
%     % x interp
%     x_r(ir  ,jr  )= 9/16*(ic) + 3/16*(ic+1) + 3/16*(ic) + 1/16*(ic+1);  % SW cell
%     x_r(ir+1,jr  )= 3/16*(ic) + 9/16*(ic+1) + 1/16*(ic) + 3/16*(ic+1);  % SE cell
%     x_r(ir  ,jr+1)= 3/16*(ic) + 1/16*(ic+1) + 9/16*(ic) + 3/16*(ic+1);  % NW cell
%     x_r(ir+1,jr+1)= 1/16*(ic) + 3/16*(ic+1) + 3/16*(ic) + 9/16*(ic+1);  % NE cell
%
%     % y interp
%     y_r(ir  ,jr  )= 9/16*(jc) + 3/16*(jc) + 3/16*(jc+1) + 1/16*(jc+1);  % SW cell
%     y_r(ir+1,jr  )= 3/16*(jc) + 9/16*(jc) + 1/16*(jc+1) + 3/16*(jc+1);  % SE cell
%     y_r(ir  ,jr+1)= 3/16*(jc) + 1/16*(jc) + 9/16*(jc+1) + 3/16*(jc+1);  % NW cell
%     y_r(ir+1,jr+1)= 1/16*(jc) + 3/16*(jc) + 3/16*(jc+1) + 9/16*(jc+1);  % NE cell
%
%     % v interp
%     v_r(ir  ,jr  )= 9/16*nc_c_h(ic,jc) + 3/16*nc_c_h(ic+1,jc) + 3/16*nc_c_h(ic,jc+1) + 1/16*nc_c_h(ic+1,jc+1);  % SW cell
%     v_r(ir+1,jr  )= 3/16*nc_c_h(ic,jc) + 9/16*nc_c_h(ic+1,jc) + 1/16*nc_c_h(ic,jc+1) + 3/16*nc_c_h(ic+1,jc+1);  % SE cell
%     v_r(ir  ,jr+1)= 3/16*nc_c_h(ic,jc) + 1/16*nc_c_h(ic+1,jc) + 9/16*nc_c_h(ic,jc+1) + 3/16*nc_c_h(ic+1,jc+1);  % NW cell
%     v_r(ir+1,jr+1)= 1/16*nc_c_h(ic,jc) + 3/16*nc_c_h(ic+1,jc) + 3/16*nc_c_h(ic,jc+1) + 9/16*nc_c_h(ic+1,jc+1);  % NE cell
%
%   end
% end
%
%% Interpolate using buffer formula
% for ic = icmin : icmax-1                   % loop as per roms indices
%   for jc = jcmin : jcmax-1
%     % Buff indices in matlab here and in roms, var indices are matlab here but just remove +1 for roms indices.
%     irm = 2 * ic +1;                         % As per matlab indices = refined grid point index for equiv. coarse even node, immediate top-right of ic
%     jrm = 2 * jc +1;
% %     ic = icl + 1;                           % As per matlab indices
% %     jc = jcl + 1;
%     shft_c = 1 + jc * Lmc;
%
%     % v interp
%     interp_r(irm  ,jrm  )= 9/16*buff(ic+shft_c) + 3/16*buff(ic+shft_c+1) + 3/16*buff(ic+shft_c+Lmc) + 1/16*buff(ic+shft_c+Lmc+1);  % SW cell
%     interp_r(irm+1,jrm  )= 3/16*buff(ic+shft_c) + 9/16*buff(ic+shft_c+1) + 1/16*buff(ic+shft_c+Lmc) + 3/16*buff(ic+shft_c+Lmc+1);  % SE cell
%     interp_r(irm  ,jrm+1)= 3/16*buff(ic+shft_c) + 1/16*buff(ic+shft_c+1) + 9/16*buff(ic+shft_c+Lmc) + 3/16*buff(ic+shft_c+Lmc+1);  % NW cell
%     interp_r(irm+1,jrm+1)= 1/16*buff(ic+shft_c) + 3/16*buff(ic+shft_c+1) + 3/16*buff(ic+shft_c+Lmc) + 9/16*buff(ic+shft_c+Lmc+1);  % NE cell
%
%   end
% end

%% Interpolate full domain
LLmc = size(nc_full_c_h,1);
MMmc = size(nc_full_c_h,2);
icmin = 0;
icmax = LLmc-1;                                   % -1 for rho numbering from 0
jcmin = 0;
jcmax = MMmc-1;                               % -1 for rho numbering from 0
buff = reshape(nc_full_c_h,1,[]); % Note no transpose.
for ic = icmin : icmax-1                   % loop as per roms indices
  for jc = jcmin : jcmax-1
    % Buff indices in matlab here and in roms, var indices are matlab here but just remove +1 for roms indices.
    irm = 2 * ic +1;                         % As per matlab indices = refined grid point index for equiv. coarse even node, immediate top-right of ic
    jrm = 2 * jc +1;
    shft_c = 1 + jc * LLmc;

    % v interp
    interp_F(irm  ,jrm  )= 9/16*buff(ic+shft_c) + 3/16*buff(ic+shft_c+1) + 3/16*buff(ic+shft_c+LLmc) + 1/16*buff(ic+shft_c+LLmc+1);  % SW cell
    interp_F(irm+1,jrm  )= 3/16*buff(ic+shft_c) + 9/16*buff(ic+shft_c+1) + 1/16*buff(ic+shft_c+LLmc) + 3/16*buff(ic+shft_c+LLmc+1);  % SE cell
    interp_F(irm  ,jrm+1)= 3/16*buff(ic+shft_c) + 1/16*buff(ic+shft_c+1) + 9/16*buff(ic+shft_c+LLmc) + 3/16*buff(ic+shft_c+LLmc+1);  % NW cell
    interp_F(irm+1,jrm+1)= 1/16*buff(ic+shft_c) + 3/16*buff(ic+shft_c+1) + 3/16*buff(ic+shft_c+LLmc) + 9/16*buff(ic+shft_c+LLmc+1);  % NE cell

  end
end
diff_F = interp_F - nc_full_r_h; % matlab interp vs raw refined
% diff_mat_roms = interp_F - nc_full_roms_int; % Difference between roms and matlab interp
% diff_roms_ref = nc_full_roms_int - nc_full_r_h; % Difference between roms and refined

%% Plot variables on grid

% figure(1)
% subplot(2,2,1)
% contourf( nc_r_h' ) ; colorbar ; title('nc-r-h') ; % Note transpose
% subplot(2,2,2)
% contourf( nc_c_h' ) ; colorbar ; title('nc-c-h') ;
% subplot(2,2,3)
% contourf( interp_r' ) ; colorbar ; title('interp-r') ;
% subplot(2,2,4)
% contourf( v_r' ) ; colorbar ; title('v-r') ;

%% Full plots
figure(2)
subplot(2,2,1)
contourf( nc_full_r_h' ) ; colorbar ; title('nc-full-r-h') ; caxis([-16 15]) ;
subplot(2,2,2)
contourf( nc_full_c_h' ) ; colorbar ; title('nc-full-c-h') ; caxis([-16 15]) ;
subplot(2,2,3)
contourf( interp_F' ) ; colorbar ; title('interp-F') ; caxis([-16 15]) ;
subplot(2,2,4)
contourf( diff_F' ) ; colorbar ; title('diff-F') ; caxis([-1 1]) ;

%% Write diff to netcdf file
file_output = '/home/devin/Documents/UCLA_2020_3/matlab/matlab-uwnd-diff-F-12-25.nc';
nccreate(file_output,'diff-F','Dimensions',{'dim1',(LLmc*2-2),'dim2',(MMmc*2-2),'dim3',Inf});
ncwrite(file_output,'diff-F',diff_F);

%% Diff plots
% figure(3)
% subplot(2,2,1)
% contourf( diff_mat_roms' ) ; colorbar ; caxis([-2e-6 1.5e-6]) ; title('diff-mat-roms') ;
% subplot(2,2,2)
% contourf( diff_roms_ref' ) ; colorbar ; caxis([-2e-3 1.5e-3]) ; title('diff-roms-ref') ;
% subplot(2,2,4)
% contourf( diff_F' ) ; colorbar ; caxis([-2e-3 1.5e-3]) ; title('diff-F') ;

% figure(4)
% subplot(2,2,1)
% contourf( diff_mat_roms' ) ; colorbar ; title('interp-F') ;
% subplot(2,2,2)
% contourf( diff_F' ) ; colorbar ; title('diff-F') ;
% subplot(2,2,3)
% contourf( interp_F' ) ; colorbar ; title('interp-F') ;
% subplot(2,2,4)
% contourf( diff_F' ) ; colorbar ; title('diff-F') ;

%% Check 1st boundary interps
% corner_c_SW = nc_c_h(1,1);
% corner_c_SE = nc_c_h(end,1);
% corner_c_NW = nc_c_h(1,end);
% corner_c_NE = nc_c_h(end,end);

%% Check sides 1st 1D interps
% side_south_r_21    = 3/4*nc_c_h(1,1)   + 1/4*nc_c_h(2,1);
% side_south_r_31    = 1/4*nc_c_h(1,1)   + 3/4*nc_c_h(2,1);
% side_west_r_12     = 3/4*nc_c_h(1,1)   + 1/4*nc_c_h(1,2);
% side_west_r_13     = 1/4*nc_c_h(1,1)   + 3/4*nc_c_h(1,2);
% side_east_r_imax2  = 3/4*nc_c_h(end,1) + 1/4*nc_c_h(end,2);
% side_east_r_imax3  = 1/4*nc_c_h(end,1) + 3/4*nc_c_h(end,2);
% side_north_r_2jmax = 3/4*nc_c_h(1,end) + 1/4*nc_c_h(2,end);
% side_north_r_3jmax = 1/4*nc_c_h(1,end) + 3/4*nc_c_h(2,end);

