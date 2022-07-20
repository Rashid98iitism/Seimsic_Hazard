%% UNCERTAINITY IN MAGNITUDE
% Inputs- 1. max_s=maximum estimated Mw for the fault s
% Output- pdf_m= PDF of uncertanity in 

max_s=load('Input/RegionB/max_mw_faults.csv');
max_s=max_s(:,1);
k=1;b=0.6136*log(10);

for ii=1:length(minnR)  
    m=4:0.1:max_s(ii); % max_s is the estimated maximum Mw of fault s 
    for i=1:length(m)
         pdf_mag(i,ii)=(b*exp(-b*(m(i)-4)))/(1-exp(-b*(max_s(ii)-4)));
    end
end

%%
b=6.136;j=1;
[~,~,data]=xlsread('regionB_EQ.xlsx');
 mag=cell2mat(data(:,1));
 
for i=1:length(mag)
    val(j)=0.02*((exp(-b*(mag(i)-4))-exp(-b*(7.1-4)))/(1-exp(-b*(7.1-4))));
    j=j+1;
end