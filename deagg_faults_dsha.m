%% DEAGGREAGATION OF SOURCES FOR DSHA 

% Inputs:- length of faults(len_s), maximum observed mag(m_obsv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BY RASHID SHAMS (11-MAY-2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of RLD

for i=1:length(m_obsv) 
     log10_RLD(i)=0.59*m_obsv(i)-2.44;
     RLD(i)=10.^(log10_RLD(i));
end


%% Mw_1: Wells and Coppersmith (1994)

for ii=1:length(log10_RLD)
    Mw_1(ii)=5.08+1.16*(log10_RLD(ii));
end

%% Mw_2: Nowroozi (1985)

  for iii=1:length(log10_RLD)
      Ms(iii)=1.259+1.244*log10_RLD(iii);
      
      % Conversion of Ms to Mw scale
      if Ms(iii)<=6.1
         Mw_2(iii)=0.67*Ms(iii)+2.12;
      else
         Mw_2(iii)=1.06*Ms(iii)-0.38;
      end
  end
  
 %% Mw_3: Regional Rupture Characterstics
 
 for k=1:length(len_s)
     pfr(k)=(RLD(k)/len_s(k))*100;
 end
    % by power law fitting
    rld_power=((121.17.*(len_s).^(-0.58)).*len_s)./100;
    
    Mw_3=(log10(rld_power)+2.44)./0.59;

%% maximum Mw from all 3 methods
for m=1:length(Mw_1)
    
    if Mw_1(m)>Mw_2(m) && Mw_1(m)>Mw_3(m)
       Mw_max(m)=Mw_1(m);
   elseif Mw_2(m)>Mw_1(m) && Mw_2(m)>Mw_3(m)
          Mw_max(m)=Mw_2(m);
   else
        Mw_max(m)=Mw_3(m);
    end
end

% Saving results in a .csv file named deagg_output_dsha.csv
output=[RLD' pfr' len_s' Mw_1' Mw_2' Mw_3' Mw_max'];
csvwrite('deagg_output.csv',output);





