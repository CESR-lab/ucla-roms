function [pf,pu,t0,phase_mkB]=egbert_correc(mjd,hr,minute,second)
%---------------------------------------------------------------------
%  Correct phases and amplitudes for real time runs
%  Use parts of pos-processing code from Egbert's & Erofeeva's (OSU) 
%  TPXO model. Their routines have been adapted from code by Richard Ray 
%  (@?) and David Cartwright.
%---------------------------------------------------------------------
%
rad=pi/180.0;
deg=180.0/pi;
tstart=mjd+hr/24+minute/(60*24)+second/(60*60*24);
disp(['Start date for nodal correction : ',mjd2greg(mjd)]);
%
% Determine nodal corrections pu & pf :
% these expressions are valid for period 1990-2010 (Cartwright 1990).
% reset time origin for astronomical arguments to 4th of May 1860:
% 
timetemp=tstart-51544.4993;
%	
% mean longitude of lunar perigee
% -------------------------------
P =  83.3535 +  0.11140353 * timetemp;
P = mod(P,360.0);
P(P<0.0) = P(P<0.0) + 360.0;
P=P*rad;
%	
% mean longitude of ascending lunar node
% --------------------------------------
N = 125.0445 -  0.05295377 * timetemp;
N = mod(N,360.0) ;
N(N<0.0) = N(N<0.0) + 360.0;
N=N*rad;
%
% nodal corrections: pf = amplitude scaling factor [], 
%                    pu = phase correction [deg]
sinn = sin(N);
cosn = cos(N);
sin2n = sin(2*N);
cos2n = cos(2*N);
sin3n = sin(3*N);
tmp1  = 1.36*cos(P)+.267*cos((P-N)); 
tmp2  = 0.64*sin(P)+.135*sin((P-N));  
temp1 = 1.-0.25*cos(2*P)-0.11*cos((2*P-N))-0.04*cosn ;
temp2 =    0.25*sin(2*P)+0.11*sin((2*P-N))+0.04*sinn ;
pftmp  = sqrt((1.-.03731*cosn+.00052*cos2n)^2+ ...
                 (.03731*sinn-.00052*sin2n)^2);% 2N2

pf( 1) = pftmp;                                                          % M2
pf( 2) = 1.0;                                                            % S2
pf( 3) = pftmp;                                                          % N2
pf( 4) = sqrt((1.+.2852*cosn+.0324*cos2n)^2+(.3108*sinn+.0324*sin2n)^2) ;% K2
pf( 5) = sqrt((1.+.1158*cosn-.0029*cos2n)^2+(.1554*sinn-.0029*sin2n)^2) ;% K1
pf( 6) = sqrt((1.+.189*cosn-0.0058*cos2n)^2+(.189*sinn -.0058*sin2n)^2) ;% O1
pf( 7) = 1.0;                                                            % P1
pf( 8) = sqrt((1.+.188*cosn)^2+(.188*sinn)^2);                           % Q1
pf( 9) = 1.043 + 0.414*cosn;                                             % Mf
pf(10) = 1.0 - 0.130*cosn ;                                              % Mm
pf(11) = pftmp^2 ;                                                       % M4
pf(12) = pftmp^2 ;                                                       % Mn4
pf(13) = pftmp^2 ;                                                       % Ms4
pf(14) = pftmp ;                                                         % 2n2
pf(15) = 1.0 ;                                                           % S1

putmp  = atan((-.03731*sinn+.00052*sin2n)/ ...
            (1.-.03731*cosn+.00052*cos2n))*deg;% 2N2

pu( 1) = putmp;                                                          % M2
pu( 2) = 0.0;                                                            % S2
pu( 3) = putmp;                                                          % N2
pu( 4) = atan(-(.3108*sinn+.0324*sin2n)/(1.+.2852*cosn+.0324*cos2n))*deg;% K2
pu( 5) = atan((-.1554*sinn+.0029*sin2n)/(1.+.1158*cosn-.0029*cos2n))*deg;% K1
pu( 6) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n;                              % O1
pu( 7) = 0.0;                                                            % P1
pu( 8) = atan(.189*sinn/(1.+.189*cosn))*deg;                             % Q1
pu( 9) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n;                             % Mf
pu(10) = 0.0;                                                            % Mm
pu(11) = putmp*2.0 ;                                                     % M4
pu(12) = putmp*2.0 ;                                                     % Mn4
pu(13) = putmp ;                                                         % Ms4
pu(14) = putmp ;                                                         % 2n2
pu(15) = 0.0 ;                                                           % S1

% to determine phase shifts below time should be in hours
% relatively Jan 1 1992 (=48622mjd) 
      
t0=48622.0*24.0;
	
% Astronomical arguments, obtained with Richard Ray's
% "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time

phase_mkB=[1.731557546,...   % M2
           0.000000000,...   % S2
           6.050721243,...   % N2
           3.487600001,...   % K2
           0.173003674,...   % K1
           1.558553872,...   % O1
           6.110181633,...   % P1
           5.877717569,...   % Q1
           1.964021610,...   % Mm
           1.756042456,...   % Mf
           3.463115091,...   % M4
           1.499093481,...   % Mn4
           1.731557546,...   % Ms4
           4.086699633,...   % 2n2
           0.000000000]*deg; % S1
              
              
