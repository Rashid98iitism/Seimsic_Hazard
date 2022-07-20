% FOCAL DEPTH ANALYSIS (TO PLOT EPICENTRAL DISTANCE VS FOCAL DEPTH)

% BY RASHID SHAMS (31-JAN-2021)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data from the declustered catalogue
[data]=csvread('declus_kishanganj.csv');
lat=data(:,2); long=data(:,1); mag=data(:,6);depth=data(:,7);

% extarct complete data as Mc=4 and exclude data where depth data is
% unavailable( or depth=0 km)
depth_comp=depth(mag>4 & depth~=0);
lat_comp=lat(mag>4 & depth~=0); long_comp=long(mag>4 & depth~=0); mag_comp=mag(mag>4 & depth~=0);

% Run for loop for computing epicentral distance ( diatance between each
% event and kishanganj(27.0982,87.9450)(lat,long)
% use distance(lat1,long1,lat2,long2) function to compute distance between
% 2 coordinates
% used deg2km function to convert obtained distance in degrees to km
% saved the distance in km in epi_dist
j=1;
for i=1:length(lat_comp)
    epi_dist(j)=deg2km(distance(lat_comp(i),long_comp(i),26.0982,87.9450));
    j=j+1;
end

% plot 
a=plot(epi_dist(mag_comp>4 & mag_comp<=5),depth_comp(mag_comp>4 & mag_comp<=5),'bo','MarkerFaceColor','b','MarkerSize',5);
hold on
b=plot(epi_dist(mag_comp>5 & mag_comp<=6),depth_comp(mag_comp>5 & mag_comp<=6),'ro','MarkerFaceColor','r','MarkerSize',5);
hold on
c=plot(epi_dist(mag_comp>6),depth_comp(mag_comp>6),'go','MarkerFaceColor','g','MarkerSize',5);

% plot specifications
set(gca, 'XAxisLocation', 'top');  % to get x axis on the top 
set(gca, 'TickDir', 'out');
set (gca,'Ydir','reverse');
ylim([0 200]);
xlim([0 550]);
yticklabels({'0','5','10','15','20','30','40','50','100','150','200'});
    
% plot labels and legend   
xlabel('Epicentral distance (km)');
ylabel('Focal Depth(km)');
title('Focal Depth Analysis');
legend([a,b,c],'4<Mw<=5','5<Mw<=6','Mw>=6','Location','southwest');