function pc = equi_tide(lon,lat,nc);

% complex equilibrium tidal potential [m]
%
%
% "M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm M4 Mn4 Ms4 2N2 S1"
%
%      amplitudes and elasticity factors for:
%
%%% Amplitude (OTPS (old version : Schwiderski, 1978)
A=[0.242334 ... % M2
   0.112743 ... % S2
   0.046397 ... % N2
   0.030684 ... % K2
   0.141565 ... % K1
   0.100661 ... % O1
   0.046848 ... % P1
   0.019273 ... % Q1
   0.042041 ... % Mf
   0.022191 ... % Mm
   0.000000 ... % M4
   0.000000 ... % Mn4
   0.000000 ... % Ms4
   0.006141 ... % 2n2
   0.000764 ];  % S1  
%%% Elasticity factor 
B=[0.693 ... % M2
   0.693 ... % S2
   0.693 ... % N2
   0.693 ... % K2
   0.736 ... % K1
   0.695 ... % O1
   0.706 ... % P1
   0.695 ... % Q1
   0.693 ... % Mf
   0.693 ... % Mm
   0.693 ... % M4
   0.693 ... % Mn4
   0.693 ... % Ms4 
   0.693 ... % 2n2 
   0.693];   % S1              
% Tide species : 2 : semidiurnal equilibrium tides
%                1 : diurnal equilibrium tides
%                0 long-period equilibrium tides
%      M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mn M4 Mn4 Ms4 2n2 s1  
ityp =[ 2  2  2  2  1  1  1  1  0  0  0   0   0   2  1]  ;
% note: for now, ispec for M4 set to 0 (ispec is only used to define forcing 
%       in atgf, and this is always  0 for M4, see OTPSnc)
%
 [nx,ny] = size(lon);

d2r = pi/180;
r2d = 180/pi;
coslat2=cos(d2r*lat).^2;      % Phase arrays
sin2lat=sin(2.*d2r*lat);
%
  cI = complex(0,1);
  pcr = zeros(nx,ny,nc);
  pc = complex(pcr,pcr);
  for ic = 1:nc
    if ityp(ic)==2                    % semidiurnal   
      p_amp = A(ic)*B(ic)*coslat2;
      p_pha = -2.*lon*d2r;
    elseif ityp(ic)==1                % diurnal
      p_amp = A(ic)*B(ic)*sin2lat;
      p_pha = -lon*d2r;
    elseif ityp(ic)==0                % long-term   
%     error: p_amp = A(ic)*B(ic)*(1-1.5*coslat2);  Error
      Pamp= A(ic)*B(ic)*(0.5-1.5*coslat2);
      p_pha = 0;
    end
    pc(:,:,ic) = p_amp.*exp(-cI*p_pha);
  end
%
%
%   disp('Get total tidal potential...')
%   Ptot=Pamp.*exp(1i*Ppha*rad) - SALamp.*exp(1i*SALpha*rad);

  return
%
