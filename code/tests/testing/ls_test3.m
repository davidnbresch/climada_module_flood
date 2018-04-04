function ls_test3()
%function to test stuff... can be deleted afterwards


%test climada_centroids_scores
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

[slope,aspect,area,TWI] = climada_centroids_scores(lon,lat,elevation,0);

[~,~,~,topoTWI] = climada_centroids_scores(lon,lat,elevation,1);

% figure
% surface(TWI)
figure
surface(topoTWI)

max(topoTWI(:))


disp('hier')

end