%% DISTANCE FOR FAULTS
% Inputs- 1.lat_s=latitude of each point in faults 
%         2. long_s=longitude of each point in faults

N=40;
 for ii=1:N
     filename =sprintf('Input/RegionB/%d.csv',ii);
     [data]=csvread(filename); 
     lat_s=data(:,1);long_s=data(:,2);
     
        for i=1:length(lat_s)
        R(i,ii)=deg2km(distance(lat_s(i),long_s(i),25.2381,87.6454));   
        end
        
        R(R==0)=Inf;
        minnR(ii)=min(R(:,ii));
 end
R(R==Inf)=0;