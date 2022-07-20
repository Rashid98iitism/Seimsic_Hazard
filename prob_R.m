%% UNCERTAINITY IN DISTANCE
% Inputs- 1. len_s=length of fault 
%         2. R= distance for fault points from a grid point 
%         3. minn= minimum distance of fault
% Output- pdf_R= pdf of unceratnity in distance

len=load('Input/RegionB/length_faults.csv');
len_s=len(:,1);
[m, n]=size(R);
pdf_R=zeros(m,length(minnR));

for ii=1:length(minnR)
    for i=1:m
        
            pdf_R(i,ii)=R(i,ii)*2/(len_s(ii)*(sqrt(R(i,ii)^(2)-minnR(ii)^(2))));
         
    end
end
 pdf_R(pdf_R==Inf)=0;pdf_R(isnan(pdf_R))=0;
       
   