%% Richardson test

df_l2 = zeros(4,6);

rname = 'Data/FMbc'
for i=0:6
    fname = [rname num2str(i)];
    
    load(fname);
    if i>0
        diff = wi-wi0;
        l2_diff = sqrt(mean( diff(:).^2));
        df_l2(1,i+1) = l2_diff;        
    end
    wi0 = wi;
end
wi1 = wi;

rname = 'Data/Cbc'
for i=0:6
    fname = [rname num2str(i)];
    
    load(fname);
    if i>0
        diff = wi-wi0;
        l2_diff = sqrt(mean( diff(:).^2));
        df_l2(2,i+1) = l2_diff;        
    end
    wi0 = wi;
    
end
diff2 = diff;
wi2 = wi;


rname = 'Data/FMnbc'
for i=0:6
    fname = [rname num2str(i)];
    
    load(fname);
    if i>0
        diff = wi-wi0;
        l2_diff = sqrt(mean(diff(:).^2));
        df_l2(3,i+1) = l2_diff;        
    end
    wi0 = wi;
end
df_l2(df_l2==0) = nan;
wi3 = wi;

rname = 'Data/Zlev'
for i=0:6
    fname = [rname num2str(i)];
    
    load(fname);
    if i>0
        diff = wi-wi0;
        l2_diff = sqrt(mean( diff(:).^2));
        df_l2(4,i+1) = l2_diff;        
    end
    wi0 = wi;
end
df_l2(df_l2==0) = nan;
wi4 = wi;
err = 0*df_l2;



load 'Data/Eksig'
Ek_sig = Ekconv;
load 'Data/Ek_Cbc'
Ek_car = Ekconv;
load 'Data/Ekzl'
Ek_zl = Ekconv;
load 'Data/Ek_FMnbc'
Ek_nbc = Ekconv;

% Richardson extrapolation to get value at inf resolution
dEk = Ek_sig(2:end)-Ek_sig(1:end-1);
Ek_lim = Ek_sig(end)+(1./3)*dEk(end);

dEk = Ek_car(2:end)-Ek_car(1:end-1);
Ek_lim2 = Ek_car(end)+(1./3)*dEk(end);

plot(abs(Ek_sig-Ek_lim),'r')
hold on
plot(abs(Ek_car-Ek_lim),'b')
plot(abs(Ek_zl-Ek_lim),'k')
plot(abs(Ek_nbc-Ek_lim),'g')
hold off
set(gca,'yscale','log')





