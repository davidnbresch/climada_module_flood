function script_FPprobSet

save_dir = 'C:\Users\Simon Rölli\Desktop\data\probSets';
num_slides = 10000;
fig = 0;
phi_hwidth = 2.5;
edges_vmax = [1 4 8 11];

%%
%load data
subS = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\subData\subS_2x3m.shp');

snapS_srtm3 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_63x92m.shp');
snapS_srtm1 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_21x30m.shp');
snapS_alti3d = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_7x10m.shp');
dem_path_srtm3 = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm3'; %srtm3
dem_path_srtm1 = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm1'; %srtm1
dem_path_alti3d = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d';%alti3d ca. 7*10m

%SRTM3
centroids = load(dem_path_srtm3,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm3 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm3 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm3 = reshape(centroids.elevation_m,n_lat,n_lon);

centroids = load(dem_path_srtm1,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm1 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm1 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm1 = reshape(centroids.elevation_m,n_lat,n_lon);

centroids = load(dem_path_alti3d,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_alti3d = reshape(centroids.lon,n_lat,n_lon);
lat_alti3d = reshape(centroids.lat,n_lat,n_lon);
elevation_alti3d = reshape(centroids.elevation_m,n_lat,n_lon);

%%
%remove some of slides in subS --> not considered in further calculations

%set some thresholds
max_area = 30000; %[m^2] according to Millage etal area (10^1–10^4 m2).
max_lgt = 1000;

%remove to large slides and others which are not covered by grid
lgt = [subS.length];
area = [subS.area];
rmv = [subS.removed];
rmv(lgt>max_lgt) = 1;
rmv(area>max_area) = 1;

%calculate probabilistic vmax_phi set
%read in angle of reach
reachAngle = [subS.reachAngle];
reachAngle(rmv~=0) = nan;

%write back rmv and reachAngle in subS
c = num2cell(rmv);
[subS.removed] = c{:};
c = num2cell(reachAngle);
[subS.reachAngle] = c{:};
%%
%define structure of phi distribution --> from histogram

[count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');

%%
%propagate slides for probabilistic parameter sets (distribution of phi not
%shifted and no cutting --> just from distribution of angle of reach
demstr = ["SRTM3" "SRTM1" "ALTI3D"];
if 0

    %srtm3
    [modS,obsS] = climada_ls_probAss(lon_srtm3,lat_srtm3,elevation_srtm3,subS,snapS_srtm3,edges_vmax,edges_phi,count_phi,num_slides,0);
    [file,path] = uiputfile('*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(1)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    %boxplot([[modS.length];[obsS.length]]')

    %srtm1
    [modS,obsS] = climada_ls_probAss(lon_srtm1,lat_srtm1,elevation_srtm1,subS,snapS_srtm1,edges_vmax,edges_phi,count_phi,num_slides,0);
    [file,path] = uiputfile('*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(2)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);

    %alti3d
    [modS,obsS] = climada_ls_probAss(lon_alti3d,lat_alti3d,elevation_alti3d,subS,snapS_alti3d,edges_vmax,edges_phi,count_phi,num_slides,0);
    [file,path] = uiputfile('*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(3)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
end

[file,path] = uigetfile('mod*.shp','load file',[save_dir filesep 'modS_' char(demstr(1)) '.shp']);
modS_SRTM3 = shaperead([path file]);
obsS_SRTM3 = shaperead([path strrep(file,'modS','obsS')]);

[file,path] = uigetfile('mod*.shp','load file',[save_dir filesep 'modS_' char(demstr(2)) '.shp']);
modS_SRTM1 = shaperead([path file]);
obsS_SRTM1 = shaperead([path strrep(file,'modS','obsS')]);

[file,path] = uigetfile('mod*.shp','load file',[save_dir filesep 'modS_' char(demstr(3)) '.shp']);
modS_ALTI3D = shaperead([path file]);
obsS_ALTI3D = shaperead([path strrep(file,'modS','obsS')]);

modS.(demstr(1)) = modS_SRTM3; obsS.(demstr(1)) = obsS_SRTM3;
modS.(demstr(2)) = modS_SRTM1; obsS.(demstr(2)) = obsS_SRTM1;
modS.(demstr(3)) = modS_ALTI3D; obsS.(demstr(3)) = obsS_ALTI3D;

figure
subplot1(2,3)


disp('hier')









end