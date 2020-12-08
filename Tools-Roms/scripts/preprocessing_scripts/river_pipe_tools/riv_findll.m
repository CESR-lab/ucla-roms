function  [i,j,rflx,mask] = findll(xr,yr,x,y,x0,y0,i,j,rflx,mask,sc0,sc1);

%%   
%%   finds the index of the lower left corner of the 
%%   enclosing box around x0,y0
%%
%%   x,y are 2d arrays of psi points
%%   i0,j0 are best guesses for the index

  if nargin < 5
    i = 1;
  end
  if nargin < 6
    j = 1;
  end
 in = i;jn=j;
 [ny,nx] = size(x);  %% psi

 count = 0;
 notfound = 1;
 idir = 1;
%while notfound
 for it = 1:10000
  if idir
    %% 'left' or 'right of the line (x,y)(j,i)->(x,y)(j,i+1)  
    if j<ny
     orient = sign((y0 - y(j,i))*(x(j+1,i) - x(j,i)) - (x0 - x(j,i))*(y(j+1,i) - y(j,i)));
    else
     orient = sign((y0 - y(j,i))*(x(j,i) - x(j-1,i)) - (x0 - x(j,i))*(y(j,i) - y(j-1,i)));
    end
    if orient<=0 
      if i < nx
        if j<ny
          new_orient = sign((y0 - y(j,i+1))*(x(j+1,i+1) - x(j,i+1)) - (x0 - x(j,i+1))*(y(j+1,i+1) - y(j,i+1)));
        else
          new_orient = sign((y0 - y(j,i+1))*(x(j,i+1) - x(j-1,i+1)) - (x0 - x(j,i+1))*(y(j,i+1) - y(j-1,i+1)));
        end
        if new_orient<=0
          i = i+1;
          count = 0;
          imax = 0;
        else
%         disp('we are in place')
          idir = 0;
          count = count+1;
          imax = 0;
        end
      else
%      disp('i is imax')
       idir = 0;
       count = count+1;
       imax = 1;
      end
    else
      imax = 0;
      if i > 1
        i = i-1;
        count = 0;
      else
%      disp('i is imin')
       idir = 0;
       count = count+1;
      end
    end
  else
    if i<nx
     orient = sign((y0 - y(j,i))*(x(j,i+1) - x(j,i)) - (x0 - x(j,i))*(y(j,i+1) - y(j,i)));
    else
     orient = sign((y0 - y(j,i))*(x(j,i) - x(j,i-1)) - (x0 - x(j,i))*(y(j,i) - y(j,i-1)));
    end
    if orient>=0 
      if j < ny
        if i<nx
          new_orient = sign((y0 - y(j+1,i))*(x(j+1,i+1) - x(j+1,i)) - (x0 - x(j+1,i))*(y(j+1,i+1) - y(j+1,i)));
        else
          new_orient = sign((y0 - y(j+1,i))*(x(j+1,i) - x(j+1,i-1)) - (x0 - x(j+1,i))*(y(j+1,i) - y(j+1,i-1)));
        end
        if new_orient >= 0
          j = j+1;
          count = 0;
          jmax = 0;
        else
%         disp('we are in place')
          idir = 1;
          count = count+1;
          jmax = 0;
        end
      else
%       disp('j is jmax')
        idir = 1;
        count = count+1;
        jmax = 1;
      end
    else
      jmax = 0;
      if j > 1
        j = j-1;
        count = 0;
      else
%      disp('j is imin')
       idir = 1;
       count = count+1;
      end
    end
  end
  if count>=2
    if imax | jmax
%    disp('outside')
     return
    end
%   disp('likely inside')
    xv = [x(j,i) x(j+1,i) x(j+1,i+1) x(j,i+1) x(j,i)];
    yv = [y(j,i) y(j+1,i) y(j+1,i+1) y(j,i+1) y(j,i)];
    if inpolygon(x0,y0,xv,yv)
%     disp('inside')
      iriv = sc1 - 1;
      hold on;plot(xr(j+1,i+1),yr(j+1,i+1),'*m');hold off
      if rflx(j+1,i+1)>0
        rflx(j+1,i+1) = 0;
        
      else
        rflx(j+1,i+1) = iriv;
      end
%     mask(j+1,i+1) = ~mask(j+1,i+1);
      nmask = 0*mask;
      nmask(~mask) = -1e5;
      nmask(rflx>0) = 0;
      hold on;pcolor(x(j:j+1,i:i+1),y(j:j+1,i:i+1),rflx(j+1:j+2,i+1:i+2)+nmask(j+1:j+2,i+1:i+2));shading flat;caxis([sc0 sc1]);hold off
    else
%     disp('outside')
    end
   return
  end
 end
 error('no convergence')



