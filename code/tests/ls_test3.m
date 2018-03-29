function ls_test3()
%function to test stuff... can be deleted afterwards

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

elevation = deminpaint(elevation);
elevation = fillsinks(elevation);

v_max = 4;


source_area = zeros(size(elevation));
source_area(100,100) = 1;

mult_flow = climada_ls_multipleflow(lon,lat,elevation);
[~,hor_dist,ver_dist] = climada_centroids_gradients(lon,lat,elevation);

spreaded = climada_ls_spread_v2(source_area,mult_flow,... 
    hor_dist,ver_dist,v_max,phi,delta_i,perWt);


disp('hier')


end