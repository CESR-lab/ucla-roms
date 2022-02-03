function create_frc_tides(gname,fname,ntides,om,constituents)

lon = ncread(gname,'lon_rho');
[nx,ny] = size(lon);


  nccreate(fname,'omega','dimensions',{'ntides',ntides},'datatype','double');
  ncwriteatt(fname,'omega','long_name','Tidal frequencies');
  ncwriteatt(fname,'omega','units','1/Second');
  ncwrite(fname,'omega',om);

  nccreate(fname,'ssh_Re','dimensions',{'xi_rho',nx,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'ssh_Re','long_name','Tidal elevation, complex; real part');
  ncwriteatt(fname,'ssh_Re','units','Meter');

  nccreate(fname,'ssh_Im','dimensions',{'xi_rho',nx,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'ssh_Im','long_name','Tidal elevation, complex; imaginary part');
  ncwriteatt(fname,'ssh_Im','units','Meter');

  nccreate(fname,'u_Re','dimensions',{'xi_u',nx-1,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'u_Re','long_name','x-direction tidal barotropic velocity, complex; real part');
  ncwriteatt(fname,'u_Re','units','Meter/Second');

  nccreate(fname,'u_Im','dimensions',{'xi_u',nx-1,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'u_Im','long_name','x-direction tidal barotropic velocity, complex; imaginary part');
  ncwriteatt(fname,'u_Im','units','Meter/Second');

  nccreate(fname,'v_Re','dimensions',{'xi_rho',nx,'eta_v',ny-1,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'v_Re','long_name','y-direction tidal barotropic velocity, complex; real part');
  ncwriteatt(fname,'v_Re','units','Meter/Second');

  nccreate(fname,'v_Im','dimensions',{'xi_rho',nx,'eta_v',ny-1,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'v_Im','long_name','y-direction tidal barotropic velocity, complex; imaginary part');
  ncwriteatt(fname,'v_Im','units','Meter/Second');

  nccreate(fname,'pot_Re','dimensions',{'xi_rho',nx,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'pot_Re','long_name','Tidal potential, complex; real part');
  ncwriteatt(fname,'pot_Re','units','Meter');

  nccreate(fname,'pot_Im','dimensions',{'xi_rho',nx,'eta_rho',ny,'ntides',ntides},'datatype','single');
  ncwriteatt(fname,'pot_Im','long_name','Tidal potential, complex; imaginary part');
  ncwriteatt(fname,'pot_Im','units','Meter');

% Constituent names of the tidal components
  ncwriteatt(fname,'/','type','Tidal elevation, velocities and potential');
  ncwriteatt(fname,'/','version','TPXO9.v2a 2020: deep=TPXO9.v1, shallow - averaged TPXO9-atlas-v2');
  ncwriteatt(fname,'/','Constituents',constituents);
