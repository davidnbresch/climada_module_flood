function ls_test3()
%function to test stuff... can be deleted afterwards

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%test ls_susceptilbitly (with ls_propagation and ls_spread(old version));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

source_area = logical(zeros(n_lat,n_lon,hazard.event_count));
for i = 1:hazard.event_count
    source_area(:,:,i) = reshape(hazard.intensity(i,:),n_lat,n_lon);
end

source = source_area(:,:,1);
exponent = 25;
dH = 1;
v_max = 4;
phi = 22;
friction = 1;
delta_i = 0.0001;
%perWt = [1 1 1 1 1 1 1 1];
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];
d2s = 1;


spread = climada_ls_susceptibility(lon,lat,elevation,source,exponent,dH,v_max,phi,friction,delta_i,perWt,d2s);

surface(elevation,spread(:,:,1));

disp('hier')

end