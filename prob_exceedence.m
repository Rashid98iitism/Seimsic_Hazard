% PROBABAILITY OF EXCEEDENCE 

x=0.02:0.01:0.5;
ns=load('Input/RegionB/Ns.csv');
N=ns(:,1);
final_ex=zeros(1,length(x));
prob_ex=zeros(length(pdf_R),length(pdf_mag));
[mm, n]=size(pdf_R);[mmm,nn]=size(pdf_mag);

for i=1:length(x)
    for ii=1:length(minnR)
        for iii=1:mm
            for iv=1:mmm
                prob_ex(iv,iii,ii)=N(ii)*pdf_mag(iv,ii)*pdf_R(iii,ii)*probb_pga(iii,iv,ii,i);
            end
        end
    end
    final_ex(i)=sum(prob_ex(:,:,:),'All');
end



%% response spectrum for hard soil site according to IS:1893 PART:3
T=0:0.01:4.0;j=1;
for i=1:length(T)
    if T(i)>=0 && T(i)<=0.10
        sa1(j)=1+15*T(i);
    elseif T(i)>=0.10 && T(i)<=0.40
        sa1(j)=2.50;
    elseif T(i)>=0.40 && T(i)<=1.0
        sa1(j)=1.0/T(i);
    end
    j=j+1;
end

%% IS CODE RESPONSE SPECTRUM for ZONE-3 (Dhanbad)
 % response spectrum for hard soil site according to IS:1893 PART:3
T=0:0.01:2.0;j=1;
for i=1:length(T)
    if T(i)>=0 && T(i)<=0.10
        sa1(j)=(1+15*T(i))*0.24;
    elseif T(i)>=0.10 && T(i)<=0.40
        sa1(j)=2.50*0.24;
    elseif T(i)>=0.40 && T(i)<=2.0
        sa1(j)=(1.0/T(i))*0.24;
    end
    j=j+1;
end
k=1;
 for i=1:length(T)
    if T(i)>=0 && T(i)<=0.10
        sa2(k)=(1+15*T(i))*0.12;
    elseif T(i)>=0.10 && T(i)<=0.40
        sa2(k)=2.50*0.12;
    elseif T(i)>=0.40 && T(i)<=2.0
        sa2(k)=(1.0/T(i))*0.12;
    end
    k=k+1;
end
hold on; plot(T,sa1,'--r',T,sa2,'--b','linewidth',2);

uhs_f_10=sgolayfilt(uhs_10,9,21);
hold on;plot(period_10,uhs_f_10,'b','linewidth',2);

hold on;
uhs_f_2=sgolayfilt(uhs_2,9,21);
hold on; plot(period_2,uhs_f_2,'r','linewidth',2);
set(gcf,'color','w');
legend('CLE Rock Site(IS-1893)-1,2002','MLE Rock Site(IS-1893)-1,2002','10% Probabilty of exceedance','2% Probability of exceedance');
