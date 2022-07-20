%%
R=length(lat);
j=1;
for i=1:length(lat)
   R(j)=deg2km(distance(lat(i),long(i),26.0982,87.9450));
    j=j+1;
end
%%
% for GMPE-1
% for magnitude Mw=

m=6;kk=1;
k=1;pga_gmpe1=length(mag);  pga_gmpe2=length(mag);
for ii=1:length(mag)
    
    y1=-1.283+0.544*m-1.0792*log10(R(ii)+exp(0.381*m))+log10(0.283);
    pga_gmpe1(k)=10^(y1);
    k=k+1;
   
    y2=0.56*m-0.0031*R(ii)-log10(R(ii)+0.0055*10*(exp(0.5*m)))+0.26+0.37;
    pga_gmpe2(kk)=(10^(y2))*0.01;
    kk=kk+1;
    
    
   
    
    
    
end







