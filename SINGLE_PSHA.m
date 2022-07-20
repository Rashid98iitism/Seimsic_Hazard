% SINGLE SOURCE PSHA CODE
%% DISTANCE
 for i=1:length(lat)
        R(i)=deg2km(distance(lat(i),long(i),23.7957,86.4304));  
 end
        
        R(R==0)=Inf;
        minnR=min(R(:));
 
R(R==Inf)=0;
%% PDF MAGNITUDE 
b=0.7338;
 m=4:0.1:max_s; % max_s is the estimated maximum Mw of fault s 
 for i=1:length(m)
      pdf_mag(i)=(b*exp(-b*(m(i)-4)))/(1-exp(-b*(max_s-4)));
 end
%% PDF DIST
m=length(R);
 for i=1:m
        
     pdf_R(i)=R(i)*2/(len_s*(sqrt(R(i)^(2)-minnR^(2))));
         
 end

 pdf_R(pdf_R==Inf)=0;pdf_R(isnan(pdf_R))=0;
 %% PROB PGA
 [mm]=length(R);
 x=0.02:0.01:0.5;
  m=4:0.1:max_s;
  for i=1:length(x)    
       for ii=1:mm
          
           for iii=1:length(m)      % LOGIC TREE OF GMPEs
       
             if R(ii)>300 
                if m(iii)<=5 && m(iii)>4  % h=focaldepth=2.5km
                a=[log((sqrt(R(ii)^2+(2.5)^2))/100) 0];
                axx=max(a);
                % pga=gmpe-3
                pga=(exp(-3.7438+1.0892*m(iii)+0.0098*(m(iii))^2-0.0046*(sqrt(R(ii)^2+(2.5)^2))-1.4817*log(sqrt(R(ii)^2+(2.5)^2))+0.0124*exp(0.9950*m(iii))+0.1249*log10(sqrt(R(ii)^2+(2.5)^2))*axx));
                probb_pga(ii,iii,i)=1-normcdf((log(x(i))-log(pga))/0.4094);
               elseif  m(iii)>5  % h=focaldepth=7.5km
               a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
               axx=max(a);
               % pga=gmpe-3
               pga=(exp(-3.7438+1.0892*m(iii)+0.0098*(m(iii))^2-0.0046*(sqrt(R(ii)^2+(7.5)^2))-1.4817*log(sqrt(R(ii)^2+(7.5)^2))+0.0124*exp(0.9950*m(iii))+0.1249*log10(sqrt(R(ii)^2+(7.5)^2))*axx));
               probb_pga(ii,iii,i)=1-normcdf((log(x(i))-log(pga))/0.4094);
                end        
             end
        
           if R(ii)>100 && R(ii)<=300
              if m(iii)<=5 && m(iii)>4
              a=[log((sqrt(R(ii)^2+(2.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii)^2+2.5^2);
              % for Fm (Magnitude Scaling Function)
              if m(iii)<=6.75
              Fm=-0.53804+0.28805*(m(iii)-6.75)-0.10164*(m(iii)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m(iii)-4.5))*log(R(ii))-0.01151*(R(ii)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0.37256;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(270/300))/log(180/300)-0.14;
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
              % pga=gmpe-1*0.30+gmpe-2*0.32+gmpe-3*0.22+gmpe-4*0.16
              pga=((10^(-1.283+0.544*m(iii)-1.072*log10((sqrt(R(ii)^2+2.5^2))+exp(0.381*m(iii)))+log10(0.283))))*0.30+((10^((0.56*m(iii)-0.0031*R(ii)-log10(sqrt(R(ii)^2+2.5^2)+0.0055*10^(0.5*m(iii)))+0.26)))/(100*9.81))*0.32+(exp(-3.7438+1.0892*m(iii)+0.0098*(m(iii))^2-0.0046*(sqrt(R(ii)^2+(2.5)^2))-1.4817*log(sqrt(R(ii)^2+(2.5)^2)+0.0124*exp(0.9950*m(iii)))+0.1249*log10(sqrt(R(ii)^2+(2.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              probb_pga(ii,iii,i)=1-normcdf((log(x(i))-log(pga))/0.3839);
              
             elseif m(iii)>5 
              a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii)^2+7.5^2);
              
              % for Fm (Magnitude Scaling Function)
              if m(iii)<=6.75
              Fm=-0.53804+0.28805*(m(iii)-6.75)-0.10164*(m(iii)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m(iii)-4.5))*log(R(ii))-0.01151*(R(ii)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0.37256;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(270/300))/log(180/300)-0.14;
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
              
              % pga=gmpe-1*0.30+gmpe-2*0.32+gmpe-3*0.22+gmpe-4*0.16
              pga=((10^(-1.283+0.544*m(iii)-1.072*log10((sqrt(R(ii)^2+7.5^2))+exp(0.381*m(iii)))+log10(0.283))))*0.30+((10^((0.56*m(iii)-0.0031*R(ii)-log10(sqrt(R(ii)^2+7.5^2)+0.0055*10^(0.5*m(iii)))+0.26)))/(100*9.81))*0.32+(exp(-3.7438+1.0892*m(iii)+0.0098*(m(iii))^2-0.0046*(sqrt(R(ii)^2+(7.5)^2))-1.4817*log(sqrt(R(ii)^2+(7.5)^2)+0.0124*exp(0.9950*m(iii)))+0.1249*log10(sqrt(R(ii)^2+(7.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              probb_pga(ii,iii,i)=1-normcdf((log(x(i))-log(pga))/0.3839);
              end
           end
        
          if R(ii)<=100 
             a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
             axx=max(a);
             % pga=gmpe-1*0.53+gmpe-2*0.12+gmpe-3*0.35
             pga=((10^(-1.283+0.544*m(iii)-1.072*log10((sqrt(R(ii)^2+7.5^2))+exp(0.381*m(iii)))+log10(0.283))))*0.53+((10^((0.56*m(iii)-0.0031*R(ii)-log10(sqrt(R(ii)^2+7.5^2)+0.0055*10^(0.5*m(iii)))+0.26)))/(100*9.81))*0.12+(exp(-3.7438+1.0892*m(iii)+0.0098*(m(iii))^2-0.0046*(sqrt(R(ii)^2+(7.5)^2))-1.4817*log(sqrt(R(ii)^2+(7.5)^2)+0.0124*exp(0.9950*m(iii)))+0.1249*log10(sqrt(R(ii)^2+(7.5)^2))*axx))*0.35;
             probb_pga(ii,iii,i)=1-normcdf((log(x(i))-log(pga))/0.3376);
          end
          
           if R(ii)==0 % no calculation for R=0 
              pga=0;
              probb_pga(ii,iii,i)=0;
           end
           
            if m(iii)==0   % no calculation for m=0
              pga=0;
              probb_pga(ii,iii,i)=0;
            end
          end
       end
  end
  %% PROB EXCEEDENCE
  [mm]=length(pdf_R);[mmm]=length(pdf_mag);

for i=1:length(x)
        for iii=1:mm
            for iv=1:mmm
                prob_ex(iv,iii)=N*pdf_mag(iv)*pdf_R(iii)*probb_pga(iii,iv,i);
            end
        end
       final_ex(i)=sum(prob_ex(:,:),'All');
 end
    

