function [z,Cs] = zlevs3(h,zeta,theta_s,theta_b,hc,N,type,scoord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function z = zlevs4(h,zeta,theta_s,theta_b,hc,N,type,scoord)
%
%  this function compute the depth of rho or w points for ROMS
%
%  On Input:
%
%    type    'r': rho point 'w': w point
%    scoord     : 'old1994' (Song, 1994),
%                 'new2006' (Sasha, 2006)
%                 'new2008' bottom stretching included (Sasha, 2008)
%    alpha,beta : optional, used for 'new2008'-type s-coordinate
%
%  On Output:
%
%    z       Depths (m) of RHO- or W-points (3D matrix).
%
%  Further Information:
%  http://www.brest.ird.fr/Roms_tools/
%
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven
%  e-mail:Pierrick.Penven@ird.fr
%
%  modified by Yusuke Uchiyama, UCLA, 2008
%  further modified by Evan Mason, UCLA, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<8
  error('Not enough input arguments')
elseif scoord=='new2008'
% disp('--- using new s-coord (2008)')
  if nargin < 9
%   disp('------ default values taken: alpha=0, beta=1')
    alpha=0; beta=1;
  end
end
%
[M, L] = size(h);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
if type=='w'
  sc = ((0:N) - N) / N;
  N = N + 1;
else
  sc=((1:N)-N-0.5) / N;
end

if (scoord=='new2008');
  % new s-coordinate allowing smooth bottom refinement
  % alpha=-1: return to pure surface s-coord; -1 < alpha <\infty: bottom refinement
  % beta:
  Cs = CSF(sc,theta_s,theta_b);
else
  % for 'old1994' and 'new2006' s-coordinate
  cff1 = 1./sinh(theta_s);
  cff2 = 0.5/tanh(0.5*theta_s);
  Cs = (1.-theta_b) * cff1 * sinh(theta_s * sc)...
      + theta_b * (cff2 * tanh(theta_s * (sc + 0.5)) - 0.5);
end
%
% Create S-coordinate system: based on model topography h(i,j),
% fast-time-averaged free-surface field and vertical coordinate
% transformation metrics compute evolving depths of of the three-
% dimensional model grid.
%
z=zeros(N,M,L);
if (scoord=='old1994')
  disp('--- using old s-coord')
  hinv=1./h;
  cff=hc*(sc-Cs);
  cff1=Cs;
  cff2=sc+1;
  for k=1:N
    z0=cff(k)+cff1(k)*h;
    z(k,:,:)=z0+zeta.*(1.+z0.*hinv);
  end
else
  if scoord=='new2006'
    disp('--- using new s-coord (2006)')
  end
  hinv=1./(h+hc);
  cff=hc*sc;
  cff1=Cs;
  for k=1:N
    z(k,:,:)=zeta+(zeta+h).*(cff(k)+cff1(k)*h).*hinv;
  end
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function Cs = CSF(sc,theta_s,theta_b)

  if (theta_s > 0.0)
    csrf=(1.-cosh(theta_s*sc))/(cosh(theta_s)-1.D0);
  else
    csrf=-sc^2;
  end
  sc1=csrf+1.0;

  if (theta_b > 0.0)
    Cs =(exp(theta_b*sc1)-1.0)/(exp(theta_b)-1.0) -1.0;
  else
    Cs = csrf;
  end

 end



