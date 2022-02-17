  function wbdcb(src,event,i,j,lon,lat,lonp,latp,rflx,mask,sc0,sc1)
     if strcmp(get(src,'SelectionType'),'normal')
        set(src,'pointer','circle')
        pointer('aim');
        cp = get(gca,'CurrentPoint');
        x0 = cp(1,1);
        y0 = cp(1,2);
        i    = evalin('base','i');
        j    = evalin('base','j');
        mask = evalin('base','mask');
        rflx = evalin('base','rflx');
        [i,j,rflx,mask] = riv_findll(lon,lat,lonp,latp,x0,y0,i,j,rflx,mask,sc0,sc1);
        assignin('base','i',i);
        assignin('base','j',j);
        assignin('base','mask',mask);
        assignin('base','rflx',rflx);

        
%       xinit = cp(1,1);yinit = cp(1,2);
%       hl = line('XData',xinit,'YData',yinit,...
%       'Marker','p','color','b');
%       set(src,'WindowButtonMotionFcn',@wbmcb)
        set(src,'WindowButtonUpFcn',@wbucb)
     end
%
%       function wbmcb(src,evnt)
%          cp = get(ah,'CurrentPoint');
%          xdat = [xinit,cp(1,1)];
%          ydat = [yinit,cp(1,2)];
%          set(hl,'XData',xdat,'YData',ydat);drawnow
%       end
%  
        function wbucb(src,evnt)
           if strcmp(get(src,'SelectionType'),'alt')
              set(src,'Pointer','arrow')
%             set(src,'WindowButtonMotionFcn','')
              set(src,'WindowButtonUpFcn','')
           else
              return
           end
        end
  end
