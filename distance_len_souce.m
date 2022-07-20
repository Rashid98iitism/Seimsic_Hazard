
for ii=1:2:no_of_col
    
     for i=1:length(length(coord))
         epi_dist(i,iii)=deg2km(distance(lat(coord(:,i)),long(coord(:,i+1)),26.0,87.94));
     end 
     
end

j=1;
for i=1:length(lat_s)
    epi_dist(j)=deg2km(distance(lat_s(i),long_s(i),26.0982,87.9450));
    j=j+1;
end