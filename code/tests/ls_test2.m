function ls_test2()
%function to test stuff... can be deleted afterwards

LCCS = climada_lc_get([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05]);

%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_srtm1_hazard.mat')

%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_srtm1_centroids.mat') 

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);
intensity = logical(zeros(n_lat,n_lon,hazard.event_count));
for i = 1:hazard.event_count
    intensity(:,:,i) = reshape(hazard.intensity(i,:),n_lat,n_lon);
end

figure
surface(LCCS.x,LCCS.y,LCCS.lc);
colormap(LCCS.cmap);

disp('hier')


end