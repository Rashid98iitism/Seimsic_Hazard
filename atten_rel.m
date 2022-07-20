 % ATTENUATION BEHAVIOUR OF VARIOUS GMPEs 
 % BY RASHID SHAMS (20-FEB-2021)
 
minnR=[10:10:500];
minnR2=[10:10:300];minnR3=[10:10:200];
max_s=7.0;x=0.02:0.01:0.5;
j=1;k=1;kk=1;jj=1;

% GMPE-3 NDMA(2010)

for ii=1:length(minnR)
a=[log((sqrt(minnR(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga3(j)=(exp(-3.7438+1.0892*max_s+0.0098*(max_s)^2-0.0046*(sqrt(minnR(ii)^2+(2.5)^2))-1.4817*log(sqrt(minnR(ii)^2+(2.5)^2)+0.0124*exp(0.9950*max_s))+0.1249*log10(sqrt(minnR(ii)^2+(2.5)^2))*axx));
j=j+1;
end

% GMPE-1 ANBU2013
for ii=1:length(minnR2)
pga1(k)=(10^(-1.283+0.544*max_s-1.072*log10((sqrt(minnR2(ii)^2+7.5^2))+exp(0.381*max_s))+log10(0.283)))/9.81;
k=k+1;
end

% GMPE-2 KANNO
for ii=1:length(minnR)
    pga2(kk)=(10^((0.56*max_s-0.0031*minnR(ii)-log10(sqrt(minnR(ii)^2+2.5^2)+0.0055*10^(0.5*max_s))+0.26)))/(100*9.81);
    kk=kk+1;
end

% GMPE-4 BOORE&ATKINSON(2008)   lnpga=Fm+Fd+Fs 
for ii=1:length(minnR3)
    
    R(ii)=sqrt(minnR3(ii)^2+2.5^2);
    
    % for Fm (Magnitude Scaling Function)
    if max_s<=6.75
       Fm=-0.53804+0.28805*(max_s-6.75)-0.10164*(max_s-6.75)^2;
    else
        Fm=-0.53804;
    end
    
    % Fd (Distance Function)
    Fd=(-0.66050+0.11970*(max_s-4.5))*log(R(ii))-0.01151*(R(ii)-1);
    
    % Fs (Site Amplification Function) Fs=Flin+Fnl  
      Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=760m/s,vref=760m/s
      % Fnl (non-linear) Fnl=Fd+Fm
      pga4nl=exp(Fm+Fd);
      
      bnl=(-0.640-(-0.14)*log(760 /300))/log(180/300)-0.14;
      del_y=bnl*log(0.09/0.06);
      del_x=log(0.09/0.03);
      c=(3*del_y-bnl*del_x)/del_x^2;
      d=-(2*del_y-bnl*del_x)/del_x^3;
      
      if pga4nl<=0.03 
          Fnl=bnl*log(0.06/0.1);
      elseif  pga4nl>0.03 && pga4nl<=0.09
          Fnl=bnl*log(0.06/0.1)+c*(log(pga4nl/0.03))^2+d*(log(pga4nl/0.03))^3;
      elseif pga4nl>0.09
          Fnl=bnl*log(pga4nl/0.1);
      end
      
      Fs=Flin+Fnl;
      
     pga4(jj)=exp(Fm+Fd+Fs);
     jj=jj+1;
end