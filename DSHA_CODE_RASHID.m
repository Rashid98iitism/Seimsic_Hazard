%% DSHA CODE 

   % Written by - RASHID SHAMS (9-MARCH-2021)
   % REFERENCE-: Sinha et al.(2019), DSHA Paper
   % Output-: FINAL MAXIMUM PGA VALUE FOR A GRID POINT 
%% 1.DISTANCE FAULTS
  % Inputs- 1.lat_s=latitude of each point in faults 
  %         2. long_s=longitude of each point in faults
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=45; 
 for ii=1:N
     filename =sprintf('Input/RegionB/%d.csv',ii);
     [data]=csvread(filename); 
     lat_s=data(:,1);long_s=data(:,2);
     
        for i=1:length(lat_s)
        R(i,ii)=deg2km(distance(lat_s(i),long_s(i),25.94392943,87.84475204));   
        end
        
        R(R==0)=Inf;
        minnR(ii)=min(R(:,ii));
 end
R(R==Inf)=0;

%% 2. Using Logic-Tree Approach (FOR DSHA CODE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs- 1. Maximum Mw for each fault 2. Minimum distance between faults
%          and grid pt. 3. x=pga values 4. M1, M2 AND M3
% Output- Gives a matrix of pga values due to all faults for a grid point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m11,m22,m33 are magnitudes, m1,m2,m3 are magnitudes multiplied with their coeff.
% R is the matrix having minimum distances with grid no(n) in columns

[sources, n]=size(R);
RR=zeros(1,n);

 for iv=1:n
     for ii=1:sources 
                                  % LOGIC TREE OF GMPEs(DSHA)
               % for e>300km
               
            if R(ii,iv)>300
                
                % FOR M1
               if m11(iv)<=5 && m11(iv)>4  % h=focaldepth=2.5km
                a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
                axx=max(a);
                % pga=gmpe-3
                pga_m1(ii,iv)=(exp(-3.7438+1.0892*m1(iv)+0.0098*(m1(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))...
                    -1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2))+0.0124*exp(0.9950*m1(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx));
               
                elseif  m11(iv)>5  % h=focaldepth=7.5km
                a=[log((sqrt(R(ii,iv)^2+(7.5)^2))/100) 0];
                axx=max(a);
               % pga=gmpe-3
               pga_m1(ii,iv)=(exp(-3.7438+1.0892*m1(iv)+0.0098*(m1(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))...
                   -1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2))+0.0124*exp(0.9950*m1(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx));
               
               end  
               
                % FOR M2
                if m22(iv)<=5 && m22(iv)>4  % h=focaldepth=2.5km
                a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
                axx=max(a);
                % pga=gmpe-3
                pga_m2(ii,iv)=(exp(-3.7438+1.0892*m2(iv)+0.0098*(m2(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))...
                    -1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2))+0.0124*exp(0.9950*m2(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx));
               
               elseif  m22(iv)>5  % h=focaldepth=7.5km
               a=[log((sqrt(R(ii,iv)^2+(7.5)^2))/100) 0];
               axx=max(a);
               % pga=gmpe-3
               pga_m2(ii,iv)=(exp(-3.7438+1.0892*m2(iv)+0.0098*(m2(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))...
                   -1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2))+0.0124*exp(0.9950*m2(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx));
                end 
                
                 % FOR M3
                if m33(iv)<=5 && m33(iv)>4  % h=focaldepth=2.5km
                a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
                axx=max(a);
                % pga=gmpe-3
                pga_m3(ii,iv)=(exp(-3.7438+1.0892*m3(iv)+0.0098*(m3(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))...
                    -1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2))+0.0124*exp(0.9950*m3(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx));
               
               elseif  m33(iv)>5  % h=focaldepth=7.5km
               a=[log((sqrt(R(ii,iv)^2+(7.5)^2))/100) 0];
               axx=max(a);
               % pga=gmpe-3
               pga_m3(ii,iv)=(exp(-3.7438+1.0892*m3(iv)+0.0098*(m3(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))...
                   -1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2))+0.0124*exp(0.9950*m3(iv))+0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx));
                end  
            end
               
             % for e>100 && e<=300
        
           if R(ii,iv)>100 && R(ii,iv)<=300
               
               % for M1
              if m11(iv)<=5 && m11(iv)>4
              a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+2.5^2);
              % for Fm (Magnitude Scaling Function)
              if m11(iv)<=6.75
              Fm=-0.53804+0.28805*(m1(iv)-6.75)-0.10164*(m1(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m1(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              %Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=760m/s,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(760/300))/log(180/300)-0.14;
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
              pga_m1(ii,iv)=((10^(-1.283+0.544*m1(iv)-1.072*log10((sqrt(R(ii,iv)^2+2.5^2))+exp(0.381*m1(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m1(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+2.5^2)+0.0055*10^(0.5*m1(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m1(iv)+0.0098*(m1(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2)+0.0124*exp(0.9950*m1(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
             elseif m11(iv)>5 
              a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+7.5^2);
              
              % for Fm (Magnitude Scaling Function)
              if m11(iv)<=6.75
              Fm=-0.53804+0.28805*(m1(iv)-6.75)-0.10164*(m1(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m1(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(760/300))/log(180/300)-0.14;
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
              pga_m1(ii,iv)=((10^(-1.283+0.544*m1(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m1(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m1(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m1(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m1(iv)+0.0098*(m1(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)+0.0124*exp(0.9950*m1(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
              end
              
              % for M2
         
              if m22(iv)<=5 && m22(iv)>4
              a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+2.5^2);
              % for Fm (Magnitude Scaling Function)
              if m22(iv)<=6.75
              Fm=-0.53804+0.28805*(m2(iv)-6.75)-0.10164*(m2(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m2(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(1500/300))/log(180/300)-0.14;
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
              pga_m2(ii,iv)=((10^(-1.283+0.544*m2(iv)-1.072*log10((sqrt(R(ii,iv)^2+2.5^2))+exp(0.381*m2(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m2(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+2.5^2)+0.0055*10^(0.5*m2(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m2(iv)+0.0098*(m2(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2)+0.0124*exp(0.9950*m2(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
             elseif m22(iv)>5 
              a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+7.5^2);
              
              % for Fm (Magnitude Scaling Function)
              if m22(iv)<=6.75
              Fm=-0.53804+0.28805*(m2(iv)-6.75)-0.10164*(m2(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m2(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(760/300))/log(180/300)-0.14;
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
              pga_m2(ii,iv)=((10^(-1.283+0.544*m2(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m2(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m2(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m2(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m2(iv)+0.0098*(m2(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)+0.0124*exp(0.9950*m2(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
              end
              
              % for M3
              if m33(iv)<=5 && m33(iv)>4
              a=[log((sqrt(R(ii,iv)^2+(2.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+2.5^2)
              % for Fm (Magnitude Scaling Function)
              if m33(iv)<=6.75
              Fm=-0.53804+0.28805*(m33(iv)-6.75)-0.10164*(m3(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m33(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(760/300))/log(180/300)-0.14;
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
              pga_m3(ii,iv)=((10^(-1.283+0.544*m3(iv)-1.072*log10((sqrt(R(ii,iv)^2+2.5^2))+exp(0.381*m3(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m3(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+2.5^2)+0.0055*10^(0.5*m3(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m3(iv)+0.0098*(m3(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(2.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(2.5)^2)+0.0124*exp(0.9950*m3(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(2.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
              
              elseif m33(iv)>5 
              a=[log((sqrt(R(ii)^2+(7.5)^2))/100) 0];
              axx=max(a);
              RR(ii)=sqrt(R(ii,iv)^2+7.5^2);
              
              % for Fm (Magnitude Scaling Function)
              if m33(iv)<=6.75
              Fm=-0.53804+0.28805*(m33(iv)-6.75)-0.10164*(m33(iv)-6.75)^2;
              else
              Fm=-0.53804;
              end
    
              % Fd (Distance Function)
              Fd=(-0.66050+0.11970*(m3(iv)-4.5))*log(R(ii,iv))-0.01151*(R(ii,iv)-1);
    
              % Fs (Site Amplification Function) Fs=Flin+Fnl  
              Flin=0;   % Flin=blin*ln(vs30/vref)   blin=-0.360,vs30=270m/c,vref=760m/s
              % Fnl (non-linear) Fnl=Fd+Fm
              pga4nl=exp(Fm+Fd);
      
              bnl=(-0.640-(-0.14)*log(760/300))/log(180/300)-0.14;
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
              pga_m3(ii,iv)=((10^(-1.283+0.544*m3(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m3(iv)))+log10(0.283))))*0.30...
                  +((10^((0.56*m3(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m3(iv)))+0.26)))/(100*9.81))*0.32...
                  +(exp(-3.7438+1.0892*m3(iv)+0.0098*(m3(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)+0.0124*exp(0.9950*m3(iv)))...
                  +0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.22+exp(Fm+Fd+Fs)*0.16;
              
              end
           end
              
           % for e<=100km
           
           % for M1,M2 and M3
           if R(ii,iv)<=100 
             a=[log((sqrt(R(ii,iv)^2+(7.5)^2))/100) 0];
             %axx=max(a);
             % pga=gmpe-1*0.53+gmpe-2*0.12+gmpe-3*0.35
             pga_m1(ii,iv)=((10^(-1.283+0.544*m1(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m1(iv)))+log10(0.283))))*0.53...
                 +((10^((0.56*m1(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m1(iv)))+0.26)))/(100*9.81))*0.12...
                 +(exp(-3.7438+1.0892*m1(iv)+0.0098*(m1(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)+0.0124*exp(0.9950*m1(iv)))...
                 +0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.35;
             pga_m2(ii,iv)=((10^(-1.283+0.544*m2(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m2(iv)))+log10(0.283))))*0.53...
                 +((10^((0.56*m2(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m2(iv)))+0.26)))/(100*9.81))*0.12...
                 +(exp(-3.7438+1.0892*m2(iv)+0.0098*(m2(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)+0.0124*exp(0.9950*m2(iv)))...
                 +0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.35;
             pga_m3(ii,iv)=((10^(-1.283+0.544*m3(iv)-1.072*log10((sqrt(R(ii,iv)^2+7.5^2))+exp(0.381*m3(iv)))+log10(0.283))))*0.53...
                 +((10^((0.56*m3(iv)-0.0031*R(ii,iv)-log10(sqrt(R(ii,iv)^2+7.5^2)+0.0055*10^(0.5*m3(iv)))+0.26)))/(100*9.81))*0.12...
                 +(exp(-3.7438+1.0892*m3(iv)+0.0098*(m3(iv))^2-0.0046*(sqrt(R(ii,iv)^2+(7.5)^2))-1.4817*log(sqrt(R(ii,iv)^2+(7.5)^2)...
                 +0.0124*exp(0.9950*m3(iv)))+0.1249*log10(sqrt(R(ii,iv)^2+(7.5)^2))*axx))*0.35;
         
           end
               
           
           if R(ii,iv)==0 % no calculation for R=0 
              pga(ii,iv)=0;
             
           end
     end
                 
 end
 

 
 %%
 % response spectrum for soft soil site according to IS:1893 PART:3
T=0:0.01:4.0;j=1;
for i=1:length(T)
    if T(i)>=0 && T(i)<=0.10
        sa(j)=1+15*T(i);
    elseif T(i)>=0.10 && T(i)<=0.67
        sa(j)=2.50;
    elseif T(i)>=0.67 && T(i)<=4.0
        sa(j)=1.67/T(i);
    end
    j=j+1;
end

    