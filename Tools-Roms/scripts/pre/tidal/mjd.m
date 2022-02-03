function [mjd]=mjd(y,m,d,h) 
%
% Get the date in Modified Julian Days
%
 if nargin<4, h=0; end;
 
 y(find(m<3))=y(find(m<3))-1;
 m(find(m<3))=m(find(m<3))+12;

 A = floor(y./100);
 B = floor(A./4);
 C = 2-A+B;
 E = floor(365.25.*(y+4716));
 F = floor(30.6001.*(m+1));
 jd= C+d+h./24+E+F-1524.5;
 mjd=jd-2400000.5;

return
