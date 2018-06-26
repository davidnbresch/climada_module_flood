function script_FPprobSet

save_dir = 'C:\Users\Simon R�lli\Desktop\data\probSets';
num_slides = 10000;
fig = 0;
phi_hwidth = 2.5;
edges_vmax = [1 4 8 11];

%%
%load data
subS = shaperead('C:\Users\Simon R�lli\Desktop\data\calibration\subData\subS_2x3m.shp');

snapS_srtm3 = shaperead('C:\Users\Simon R�lli\Desktop\data\calibration\snapData\snapS_63x92m.shp');
snapS_srtm1 = shaperead('C:\Users\Simon R�lli\Desktop\data\calibration\snapData\snapS_21x30m.shp');
snapS_alti3d = shaperead('C:\Users\Simon R�lli\Desktop\data\calibration\snapData\snapS_7x10m.shp');
dem_path_srtm3 = 'C:\Users\Simon R�lli\Desktop\data\centroids_large\_LS_Sarnen_srtm3'; %srtm3
dem_path_srtm1 = 'C:\Users\Simon R�lli\Desktop\data\centroids_large\_LS_Sarnen_srtm1'; %srtm1
dem_path_alti3d = 'C:\Users\Simon R�lli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d';%alti3d ca. 7*10m

%SRTM3
centroids = load(dem_path_srtm3,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm3 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm3 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm3 = reshape(centroids.elevation_m,n_lat,n_lon);

%SRTM1
centroids = load(dem_path_srtm1,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm1 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm1 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm1 = reshape(centroids.elevation_m,n_lat,n_lon);

%ALTI3D
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
max_area = 30000; %[m^2] according to Millage etal area (10^1�10^4 m2).
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

%remove them completely
rmv = find(rmv~=0);
subS(rmv) = [];
snapS_srtm3(rmv) = [];
snapS_srtm1(rmv) = [];
snapS_alti3d(rmv) = [];
%%
%define structure of phi distribution --> from histogram

[count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');

%%
%propagate slides for probabilistic parameter sets (distribution of phi not
%shifted and no cutting --> just from distribution of angle of reach
demstr = ["SRTM3" "SRTM1" "ALTI3D"];
if 0

    %srtm3
    [modS,obsS] = climada_ls_probAss(lon_srtm3,lat_srtm3,elevation_srtm3,subS,snapS_srtm3,edges_vmax,edges_phi,count_phi,num_slides,fig);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(1)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    %boxplot([[modS.length];[obsS.length]]')

    %srtm1
    [modS,obsS] = climada_ls_probAss(lon_srtm1,lat_srtm1,elevation_srtm1,subS,snapS_srtm1,edges_vmax,edges_phi,count_phi,num_slides,0);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(2)) '.shp']);     

    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);

    %alti3d
    [modS,obsS] = climada_ls_probAss(lon_alti3d,lat_alti3d,elevation_alti3d,subS,snapS_alti3d,edges_vmax,edges_phi,count_phi,num_slides,0);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(3)) '.shp']);
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

%plots
climada_ls_probAssPlot(modS,obsS,'length','both',[0 2500],50,[0 2500])
climada_ls_probAssPlot(modS,obsS,'length','qq',[0 500])
climada_ls_probAssPlot(modS,obsS,'length','hist','',100,[0 2500])
climada_ls_probAssPlot(modS,obsS,'length','hist','',25,[0 1000])

%histogram of distribution of phi for slides lenght = 0 and slides
%vs distribtuion of vmax phi for slides >0 --> maybe you see that it is
%just dependend on phi (or this is the case by definition to be clear)

%find zero_slides
zero_idx = find([modS_SRTM3.length]==0);
modS_zeros.SRTM3 = modS_SRTM3(zero_idx);
obsS_zeros.SRTM3 = obsS_SRTM3(zero_idx);
zero_idx = find([modS_SRTM1.length]==0);
modS_zeros.SRTM1 = modS_SRTM1(zero_idx);
obsS_zeros.SRTM1 = obsS_SRTM1(zero_idx);
zero_idx = find([modS_ALTI3D.length]==0);
modS_zeros.ALTI3D = modS_ALTI3D(zero_idx);
obsS_zeros.ALTI3D = obsS_ALTI3D(zero_idx);

%find slides bigger than zero
nonzero_idx = find([modS_SRTM3.length]>0);
modS_nonzeros.SRTM3 = modS_SRTM3(nonzero_idx);
obsS_nonzeros.SRTM3 = obsS_SRTM3(nonzero_idx);
nonzero_idx = find([modS_SRTM1.length]>0);
modS_nonzeros.SRTM1 = modS_SRTM1(nonzero_idx);
obsS_nonzeros.SRTM1 = obsS_SRTM1(nonzero_idx);
nonzero_idx = find([modS_ALTI3D.length]>0);
modS_nonzeros.ALTI3D = modS_ALTI3D(nonzero_idx);
obsS_nonzeros.ALTI3D = obsS_ALTI3D(nonzero_idx);

%dist zeros vs nonzeros phi
climada_ls_probAssPlot(modS_zeros,modS_nonzeros,'phi','hist','',2.5,[0 60])
%dist zeros vs nonzeros vmax
%climada_ls_probAssPlot(modS_zeros,modS_nonzeros,'vmax','hist','',1,[1 11])


%histogram of distribution of phi for slides lenght > threshold and slides
%vs distribtuion of vmax phi for slides >0 --> maybe you see that it is
%just dependend on phi (or this is the case by definition to be clear)

str_lgt_th = '1000'; %threshold for length
lgt_th = str2num(str_lgt_th);

%find slides smaller than threshold
idx = find([modS_SRTM3.length]<=lgt_th);
modS_stLgtTh.SRTM3 = modS_SRTM3(idx);
obsS_stLgtTh.SRTM3 = obsS_SRTM3(idx);
idx = find([modS_SRTM1.length]<=lgt_th);
modS_stLgtTh.SRTM1 = modS_SRTM1(idx);
obsS_stLgtTh.SRTM1 = obsS_SRTM1(idx);
idx = find([modS_ALTI3D.length]<=lgt_th);
modS_stLgtTh.ALTI3D = modS_ALTI3D(idx);
obsS_stLgtTh.ALTI3D = obsS_ALTI3D(idx);

%find slides bigger than threshold
idx = find([modS_SRTM3.length]>lgt_th);
modS_gtLgtTh.SRTM3 = modS_SRTM3(idx);
obsS_gtLgtTh.SRTM3 = obsS_SRTM3(idx);
idx = find([modS_SRTM1.length]>lgt_th);
modS_gtLgtTh.SRTM1 = modS_SRTM1(idx);
obsS_gtLgtTh.SRTM1 = obsS_SRTM1(idx);
idx = find([modS_ALTI3D.length]>lgt_th);
modS_gtLgtTh.ALTI3D = modS_ALTI3D(idx);
obsS_gtLgtTh.ALTI3D = obsS_ALTI3D(idx);

%dist st threshold vs gt threshold for phi
climada_ls_probAssPlot(modS_stLgtTh,modS_gtLgtTh,'phi','hist','',2.5,[0 60])
%dist st threshold vs gt threshold for vmax
climada_ls_probAssPlot(modS_stLgtTh,modS_gtLgtTh,'vmax','hist','',1,[1 11])

%dist and qq when zero slides are not considered
climada_ls_probAssPlot(modS_nonzeros,obsS_nonzeros,'length','both',[0 2500],50,[0 2500])





disp('hier')









end