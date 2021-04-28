clc; clear

%% info

%  This tool enables the use of online child boundary input file from
%  parent grid simulation in roms.
%
%  This tool generates the child boundary positions for the child boundary rho,
%  u and v points in terms of the parent grid coordinates.
%  Outputs:
%  - child bry rho-points in terms of parent rho-points = ic_rho, jc_rho
%  - child bry   u-points in terms of parent   u-points = ic_u,   jc_u
%  - child bry   u-points in terms of parent   v-points = ic_v,   jc_v
%
%  On a roms grid a rh-var will have i and j positions (xi_rho, eta_rho),
%  a u-var will have (xi_u, eta_rho) and a v-var will have (xi_rho, eta_v).
%
%  This means in thr same grid a u-point will have the same j value as a 
%  rho-point if positioned on thesame row. 
%  However, since the child grid is likely not orientated with its x & y axes the
%  same way as the parent grid, the child rho-var j value in terms of the parent
%  j coordinate will not match the j value of a child u-var in term of the
%  parent coordinate!
%  Therefore, we cannot share directly i-rho and j-rho for u/v points compared 
%  with rho-points. E.g. jc_rho /= jc_u and ic_rho /= ic_v.
%
%  For i_u we +0.5 because we first average to u point from rho:
%  r0  u1  r1 -> thus u1 is at 0.5 rho, so +0.5 to get u=1.
%  The same is done for j_v.


%% init

% START USER INPUT ----------
grids_dir  = '/home/devin/code-roms-versions/eclipse/WEC/WEC_simus/bry_child_2021_03/'; % directory of grid files
pname      = append(grids_dir, 'sample_grd_riv.nc');                        % parent grid
pgrid_info = 'Sample domain - Santa Barbara. 120m grid.';                   % identifiable info about grid to add to output file

cname      = append(grids_dir, 'sample_child_grd.nc');                       % child grid
cgrid_info = 'Sample child grid - Santa Barbara. 50m grid.';

bcname     = append(grids_dir, 'sample_child_ij_bry.nc');                   % output file name of i/j child points in parent coords
delete(bcname);                                                             % prevents overwriting error

OBC_SWNE = [1,1,1,0];                                                       % choose open boundaries: south,west,north,east. 1=yes, 0=no.
% END USER INPUT ------------

bry_names  = [ "south"   "west"     "north"   "east"    ];
dim_names  = [ "xi_rho"  "eta_rho"  "xi_rho"  "eta_rho" ];
dim_unames = [ "xi_u"    "eta_rho"  "xi_u"    "eta_rho" ];
dim_vnames = [ "xi_rho"  "eta_v"    "xi_rho"  "eta_v"   ];

lonc     = ncread(cname,'lon_rho') - 360;  
latc     = ncread(cname,'lat_rho');
mskc     = ncread(cname,'mask_rho');
anglec_r = ncread(cname,'angle');

lonp = ncread(pname,'lon_rho');        
latp = ncread(pname,'lat_rho');

[nxc,nyc] = size(lonc);                                                     % dimension sizes of child domain
[nxp,nyp] = size(lonp);                                                     % dimension sizes of parent domain

b_rdimsize = [ nxc,   nyc,   nxc,   nyc   ];                                % array of child boundary rho dimension sizes for SWNE
b_udimsize = [ nxc-1, nyc,   nxc-1, nyc   ];                                % array of child boundary u-point dimension sizes for SWNE
b_vdimsize = [ nxc,   nyc-1, nxc,   nyc-1 ];                                % array of child boundary v-point dimension sizes for SWNE

b_xrows = [ 1, nxc ; 1, 2   ; 1    , nxc ; nxc-1, nxc ];                    % rows closest to each boundary
b_yrows = [ 1, 2   ; 1, nyc ; nyc-1, nyc ; 1    , nyc ];                    % but row can be i or j
b_indx  = [ 1      ; 1      ; 2          ; 2          ];                    % which of 2 rho rows to use (outer edge of boundary)

ip_rho = [0:nxp-1];                                                         % rho numbering of parent domain
jp_rho = [0:nyp-1];
[ip_rho,jp_rho] = meshgrid(ip_rho,jp_rho);

ip_rho = ip_rho';
jp_rho = jp_rho';

for b=1:4

	if OBC_SWNE(b)                                                          % if open boundary
    
        %% rho positions:
        
        irbry_name = strcat('i_rho_',bry_names(b));       
        jrbry_name = strcat('j_rho_',bry_names(b));
        anglb_name = strcat('angle_rho_',bry_names(b));

        blname = strcat(bry_names(b),' child boundary in parent grid coordinates'); % long name

        lonbc = lonc( b_xrows(b,1):b_xrows(b,2) , b_yrows(b,1):b_yrows(b,2) );      % 1st 2 rows/cols (S & W) or last 2 rows/cols (N & E) of grid
        latbc = latc( b_xrows(b,1):b_xrows(b,2) , b_yrows(b,1):b_yrows(b,2) );
        mskbc = mskc( b_xrows(b,1):b_xrows(b,2) , b_yrows(b,1):b_yrows(b,2) );
        angcr2 = anglec_r( b_xrows(b,1):b_xrows(b,2) , b_yrows(b,1):b_yrows(b,2) );

        ic_rho = griddata(lonp,latp,ip_rho,lonbc,latbc).*mskbc+(mskbc-1);           % +(mskbc-1) does nothing if ocean, and sets value to -1 if mask
        jc_rho = griddata(lonp,latp,jp_rho,lonbc,latbc).*mskbc+(mskbc-1);           
                
        ic_rho(isnan(ic_rho)) = -2;                                                 % set to -2 if not in the domain
        jc_rho(isnan(jc_rho)) = -2;
        
      
        nccreate(bcname, irbry_name, 'dimensions', {dim_names(b), b_rdimsize(b)}, 'datatype', 'single');
        ncwriteatt(bcname, irbry_name, 'long_name', strcat('i-rho coordinate of',{' '}, blname));
        ncwriteatt(bcname, irbry_name, 'units', 'non-dimensional');

        nccreate(bcname, jrbry_name, 'dimensions', {dim_names(b), b_rdimsize(b)},'datatype','single');
        ncwriteatt(bcname, jrbry_name, 'long_name', strcat('j-rho coordinate of',{' '}, blname));
        ncwriteatt(bcname, jrbry_name, 'units', 'non-dimensional');       

        if b==1 || b==3                                                             % transpose indices if W or E boundaries (since only printing one row/column)
            ncwrite(bcname, irbry_name, ic_rho( :,b_indx(b) ) );  
            ncwrite(bcname, jrbry_name, jc_rho( :,b_indx(b) ) );
        else
            ncwrite(bcname, irbry_name, ic_rho( b_indx(b),: )');
            ncwrite(bcname, jrbry_name, jc_rho( b_indx(b),: )');
        end        
        
        %% u positions:
        
        if b==1 || b==3                                                     % S/N boundaries
          
            for i=1:nxc-1                                                                     % u-points one less than rho-points in x hence -1
                ic_u(i) = 0.5 + 0.5 * ( ic_rho( i, b_indx(b) ) + ic_rho( i+1, b_indx(b) ) );  % b_indx(b) is to choose outer edge of 2 bry rows. 
                                                                                              % +0.5 for u-point shift from rho point.
                jc_u(i) =       0.5 * ( jc_rho( i, b_indx(b) ) + jc_rho( i+1, b_indx(b) ) );  % no +0.5 as j is at rho coordinate for u-points
                                                                                              % jc_u not the same as jc_rho since child is at an angle to parent.                                                                                             
                if ic_rho( i, b_indx(b) ) == -1 || ic_rho( i+1, b_indx(b) ) == -1             % catch umask! rmask was set to -1 above in rho section
                    ic_u(i) = -1;
                    jc_u(i) = -1;
                end   
                
                anglec_u(i) =   0.5 * ( angcr2( i, b_indx(b) ) + angcr2( i+1, b_indx(b) ) );  % child angle at bry u-points                
            end
        else             % W/E boundaries (transpose)                     
          
            for j=1:nyc                                                                       % u-points same number as rho-points in y
                
                ic_u(j) = 0.5 + 0.5 * ( ic_rho( 1 , j )        + ic_rho( 2, j ) );            
                jc_u(j) =       0.5 * ( jc_rho( 1 , j )        + jc_rho( 2, j ) );            % no +0.5 as j is at rho coordinate for u points
                
                if ic_rho( 1 , j ) == -1 || ic_rho( 2, j ) == -1                              % catch umask
                    ic_u(j) = -1;
                    jc_u(j) = -1;
                end
                
                anglec_u(j) =    0.5 * ( angcr2( 1 , j )       + angcr2( 2 , j ) );              
            end            
        end      
              
        iubry_name = strcat('i_u_',bry_names(b));        
        nccreate(bcname, iubry_name, 'dimensions', {dim_unames(b), b_udimsize(b)}, 'datatype', 'single');
        ncwriteatt(bcname, iubry_name, 'long_name', strcat('i-u coordinate of',{' '}, blname));
        ncwriteatt(bcname, iubry_name, 'units', 'non-dimensional');
        ncwrite(bcname, iubry_name, ic_u );
        clear ic_u;                                                         % need to delete variable since assigned in loop so size won't be correct.  
        
        jubry_name = strcat('j_u_',bry_names(b));        
        nccreate(bcname, jubry_name, 'dimensions', {dim_unames(b), b_udimsize(b)}, 'datatype', 'single');
        ncwriteatt(bcname, jubry_name, 'long_name', strcat('j-rho for u point of',{' '}, blname));
        ncwriteatt(bcname, jubry_name, 'units', 'non-dimensional');
        ncwrite(bcname, jubry_name, jc_u );
        clear jc_u;

        anglb_name = strcat('angle_u_',bry_names(b));
        nccreate(bcname, anglb_name, 'dimensions', {dim_unames(b), b_udimsize(b)},'datatype','single');
        ncwriteatt(bcname, anglb_name, 'long_name', strcat('angle at u-point of',{' '}, blname));
        ncwriteatt(bcname, anglb_name, 'units', 'radians'); 
        ncwrite(bcname, anglb_name, anglec_u );
        clear anglec_u;
                    
        %% v positions:
      
        if b==1 || b==3                                                     % S/N boundaries
          
            for i=1:nxc                                                     % v-points same number as rho-points in x
                ic_v(i) =       0.5 * ( ic_rho( i, 1 ) + ic_rho( i, 2 ) );  % no +0.5 as i is at rho coordinate for u-points
                                                                            % ic_v not the same as ic_rho since child axis is at an angle to parent axis.
                jc_v(i) = 0.5 + 0.5 * ( jc_rho( i, 1 ) + jc_rho( i, 2 ) );  % +0.5 for v-point shift in j
                                                                            
                
                if ic_rho( i, 1 ) == -1 || ic_rho( i, 2 ) == -1             % catch vmask
                    ic_v(i) = -1;
                    jc_v(i) = -1;
                end
                
                anglec_v(i) =       0.5 * ( angcr2( i, 1 ) + angcr2( i, 2 ) ); % child angle at bry v-points                
            end
        else                                                                % W/E boundaries (transpose)                     
          
            for j=1:nyc-1                                                   % v-points one less than rho-points in y hence -1
                
                ic_v(j) =       0.5 * ( ic_rho( b_indx(b) , j ) + ic_rho( b_indx(b), j+1 ) );            
                jc_v(j) = 0.5 + 0.5 * ( jc_rho( b_indx(b) , j ) + jc_rho( b_indx(b), j+1 ) ); % no +0.5 as j is at rho coordinate for u points
                
                
                if ic_rho( b_indx(b) , j ) == -1 || ic_rho( b_indx(b), j+1 ) == -1            % catch vmask
                    ic_v(j) = -1;
                    jc_v(j) = -1;
                end
                
                anglec_v(j) =   0.5 * ( angcr2( b_indx(b) , j ) + angcr2( b_indx(b), j+1 ) );
            end            
        end      
      
        ivbry_name = strcat('i_v_',bry_names(b));        
        nccreate(bcname, ivbry_name, 'dimensions', {dim_vnames(b), b_vdimsize(b)}, 'datatype', 'single');
        ncwriteatt(bcname, ivbry_name, 'long_name', strcat('i-rho for v point of',{' '}, blname));
        ncwriteatt(bcname, ivbry_name, 'units', 'non-dimensional');
        ncwrite(bcname, ivbry_name, ic_v );
        clear ic_v;
        
        jvbry_name = strcat('j_v_',bry_names(b));        
        nccreate(bcname, jvbry_name, 'dimensions', {dim_vnames(b), b_vdimsize(b)}, 'datatype', 'single');
        ncwriteatt(bcname, jvbry_name, 'long_name', strcat('j-v coordinate of',{' '}, blname));
        ncwriteatt(bcname, jvbry_name, 'units', 'non-dimensional');
        ncwrite(bcname, jvbry_name, jc_v );
        clear jc_v;     

        anglb_name = strcat('angle_v_',bry_names(b));
        nccreate(bcname, anglb_name, 'dimensions', {dim_vnames(b), b_vdimsize(b)},'datatype','single');
        ncwriteatt(bcname, anglb_name, 'long_name', strcat('angle at v-point of',{' '}, blname));
        ncwriteatt(bcname, anglb_name, 'units', 'radians'); 
        ncwrite(bcname, anglb_name, anglec_v );
        clear anglec_v;        

	end
   
end
      
ncwriteatt(bcname, '/', 'parent_grid_info', pgrid_info);                    % add global attribute to understand file grids  
ncwriteatt(bcname, '/', 'child_grid_info',  cgrid_info);
