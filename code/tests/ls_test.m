function ls_test()
%function to test stuff... can be deleted afterwards

% load('C:\Users\Simon Rölli\Desktop\data\z.mat','z');
% load('C:\Users\Simon Rölli\Desktop\data\lon.mat','lon');
% load('C:\Users\Simon Rölli\Desktop\data\lat.mat','lat');

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')

%get gridded datasets
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
z = reshape(centroids.elevation_m,n_lat,n_lon);

z = deminpaint(z);
z = fillsinks(z);

nan = isnan(z);
if any(nan(:))
    cprintf([1 0.5 0], 'WARNING: your elevation grid includes some NaNs. The calucaltion of the flow accumulation may encounter\n')
    cprintf([1 0.5 0], '\t some problems. Check the results and consider removing the NaNs (e.g. by using deminpaint(z)).')
end

%with topotoolbox
DEM = GRIDobj(lon,lat,flipud(z));
flow_acc1 = climada_ls_flowacc(DEM,1);
figure
surface(lon,lat,z,log(flow_acc1),'LineStyle','none')
colorbar
caxis([0 10])

%with climada code
exponent = 1.1;
dH = 0;
flat_areas = 1;

mult_flow = climada_ls_multipleflow(lon,lat,z,exponent,dH,flat_areas);
flow_acc2 = climada_ls_flowacc(mult_flow,0);
figure
surface(lon,lat,z,log(flow_acc2),'LineStyle','none')
colorbar
caxis([0 10])

disp('hier')


end