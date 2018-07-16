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

% figure
% boxplot([[subS.reachAngle];[snapS_alti3d.reachAngle];[snapS_srtm1.reachAngle];[snapS_srtm3.reachAngle]]',...
%     {'ALTI3D (2m)','ALTI3D (10m)','SRTM1 (30m)','SRTM3 (90m)'},'Widths',[0.2 0.2 0.2 0.2],...
%     'FactorGap',0.001)
% ylabel('Angle of Reach [\circ]')
% yl = ylim;
% ylim([0 yl(2)])
% xl = xlim; yl = ylim;
% hold on
% plot([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'-','LineWidth',2,'Color','black')

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
max_area = 30000; %[m^2] according to Millage etal area (10^1–10^4 m2).
max_lgt = 1000;

%remove to large slides and others which are not covered by grid
lgt = [subS.length];
area_dat = [subS.area];
rmv = [subS.removed];
rmv(lgt>max_lgt) = 1;
rmv(area_dat>max_area) = 1;

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
% caculation with different energy decay values    
%     path = [save_dir];
%     for en_decay = 0.1:0.1:1
%         disp(['energy decay' num2str(en_decay)])
%         %srtm3
%         [modS,obsS] = climada_ls_probAss(lon_srtm3,lat_srtm3,elevation_srtm3,subS,snapS_srtm3,edges_vmax,edges_phi,count_phi,num_slides,en_decay,fig);
%         file = ['modS_SRTM3_decay' erase(num2str(en_decay),'.') '.shp'];
%         shapewrite(modS,[path filesep file]);
%         shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
%         %srtm1
%         [modS,obsS] = climada_ls_probAss(lon_srtm1,lat_srtm1,elevation_srtm1,subS,snapS_srtm1,edges_vmax,edges_phi,count_phi,num_slides,en_decay,0);
%         file = ['modS_SRTM1_decay' erase(num2str(en_decay),'.') '.shp'];
%         shapewrite(modS,[path filesep file]);
%         shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
%         %alti3d
%         [modS,obsS] = climada_ls_probAss(lon_alti3d,lat_alti3d,elevation_alti3d,subS,snapS_alti3d,edges_vmax,edges_phi,count_phi,num_slides,en_decay,0);
%         file = ['modS_ALTI3D_decay' erase(num2str(en_decay),'.') '.shp'];
%         shapewrite(modS,[path filesep file]);
%         shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]); 
%     end
    
    en_decay = 0.9;
    %srtm3
    [modS,obsS] = climada_ls_probAss(lon_srtm3,lat_srtm3,elevation_srtm3,subS,snapS_srtm3,edges_vmax,edges_phi,count_phi,num_slides,en_decay,fig);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(1)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    %boxplot([[modS.length];[obsS.length]]')

    %srtm1
    [modS,obsS] = climada_ls_probAss(lon_srtm1,lat_srtm1,elevation_srtm1,subS,snapS_srtm1,edges_vmax,edges_phi,count_phi,num_slides,en_decay,0);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(2)) '.shp']);     

    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);

    %alti3d
    [modS,obsS] = climada_ls_probAss(lon_alti3d,lat_alti3d,elevation_alti3d,subS,snapS_alti3d,edges_vmax,edges_phi,count_phi,num_slides,en_decay,0);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_' char(demstr(3)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    
    
end

%%
%analysis when considering distribution of angle of reach for corresponding
%snap
if 0
    
    %srtm3
    fig = 1;
    en_decay = 0;
    reachAngle = [snapS_srtm3.reachAngle];
    %remove angle of reach smaller than zero 
    reachAngle(reachAngle<=0) = nan;
    [count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');
    [modS,obsS] = climada_ls_probAss(lon_srtm3,lat_srtm3,elevation_srtm3,subS,snapS_srtm3,edges_vmax,edges_phi,count_phi,num_slides,en_decay,fig);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_rA' char(demstr(1)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    
    %srtm1
    fig = 1;
    en_decay = 0;
    reachAngle = [snapS_srtm1.reachAngle];
    reachAngle(reachAngle<=0) = nan;
    [count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');
    [modS,obsS] = climada_ls_probAss(lon_srtm1,lat_srtm1,elevation_srtm1,subS,snapS_srtm1,edges_vmax,edges_phi,count_phi,num_slides,en_decay,fig);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_rA' char(demstr(2)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
    
    %alti3d
    fig = 1; 
    en_decay = 0;
    reachAngle = [snapS_alti3d.reachAngle];
    reachAngle(reachAngle<0) = nan;
    [count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');
    [modS,obsS] = climada_ls_probAss(lon_alti3d,lat_alti3d,elevation_alti3d,subS,snapS_alti3d,edges_vmax,edges_phi,count_phi,num_slides,en_decay,fig);
    [file,path] = uiputfile('mod*.shp','Save file as',[save_dir filesep 'modS_rA' char(demstr(3)) '.shp']);
    shapewrite(modS,[path filesep file]);
    shapewrite(obsS,[path filesep strrep(file,'modS','obsS')]);
end

path = 'C:\Users\Simon Rölli\Desktop\data\probSets\';

%normal data with angle of reach from 2m
modS_SRTM3 = shaperead([path 'modS_SRTM3.shp']);
obsS_SRTM3 = shaperead([path 'obsS_SRTM3.shp']);
modS_SRTM1 = shaperead([path 'modS_SRTM1.shp']);
obsS_SRTM1 = shaperead([path 'obsS_SRTM1.shp']);
modS_ALTI3D = shaperead([path 'modS_ALTI3D.shp']);
obsS_ALTI3D = shaperead([path 'obsS_ALTI3D.shp']);
%angle of reach from corresponding dem
modS_SRTM3_rA = shaperead([path 'modS_SRTM3_rA.shp']);
obsS_SRTM3_rA = shaperead([path 'obsS_SRTM3_rA.shp']);
modS_SRTM1_rA = shaperead([path 'modS_SRTM1_rA.shp']);
obsS_SRTM1_rA = shaperead([path 'obsS_SRTM1_rA.shp']);
modS_ALTI3D_rA = shaperead([path 'modS_ALTI3D_rA.shp']);
obsS_ALTI3D_rA = shaperead([path 'obsS_ALTI3D_rA.shp']);

%construct structure for plot
modS_rA.SRTM3.SRTM3 = modS_SRTM3;
modS_rA.SRTM3.SRTM3_angle = modS_SRTM3_rA;
obsS_rA.SRTM3.SRTM3 = obsS_SRTM3;
obsS_rA.SRTM3.SRTM3_angle = obsS_SRTM3_rA;
modS_rA.SRTM1.SRTM1 = modS_SRTM1;
modS_rA.SRTM1.SRTM1_angle = modS_SRTM1_rA;
obsS_rA.SRTM1.SRTM1 = obsS_SRTM1;
obsS_rA.SRTM1.SRTM1_angle = obsS_SRTM1_rA;
modS_rA.ALTI3D.ALTI3D = modS_ALTI3D;
modS_rA.ALTI3D.ALTI3D_angle = modS_ALTI3D_rA;
obsS_rA.ALTI3D.ALTI3D = obsS_ALTI3D;
obsS_rA.ALTI3D.ALTI3D_angle = obsS_ALTI3D_rA;

field = 'length'
climada_ls_probAssPlot(modS_rA,obsS_rA,field,'both',[0 2500],50,[0 2500])

%plot cdf of with modelled length with normal angle of reach vs cdf of
%modelled length with corresponding angle of reach
modS_cdf.SRTM3.SRTM3 = modS_SRTM3;
modS_cdf.SRTM3.SRTM3_rA = modS_SRTM3_rA;
modS_cdf.SRTM1.SRTM1 = modS_SRTM1;
modS_cdf.SRTM1.SRTM1_rA = modS_SRTM1_rA;
modS_cdf.ALTI3D.ALTI3D = modS_ALTI3D;
modS_cdf.ALTI3D.ALTI3D_rA = modS_ALTI3D_rA;

climada_ls_probAssPlot(modS_cdf,modS_cdf,field,'both',[0 2500],50,[0 2500],1,[0 1000])


%%

%read in files for energy decay
if 0
    [file,path] = uigetfile('mod*SRTM3*.shp','load file',[save_dir filesep 'modS_' char(demstr(1)) '.shp'],'MultiSelect','on');
    if isstr(file), file={file}; end
    for i=1:numel(file)
        struc_name = char(file(i)); struc_name = struc_name(1:end-4);
        modS_temp = shaperead([path char(file(i))]);
        modS.SRTM3.(struc_name) = modS_temp;
        obsS_temp = shaperead([path strrep(char(file(i)),'modS','obsS')]);
        obsS.SRTM3.(struc_name) = obsS_temp;
    end

    [file,path] = uigetfile('mod*SRTM1*.shp','load file',[save_dir filesep 'modS_' char(demstr(2)) '.shp'],'MultiSelect','on');
    if isstr(file), file={file}; end
    for i=1:numel(file)
        struc_name = char(file(i)); struc_name = struc_name(1:end-4);
        modS_temp = shaperead([path char(file(i))]);
        modS.SRTM1.(struc_name) = modS_temp;
        obsS_temp = shaperead([path strrep(char(file(i)),'modS','obsS')]);
        obsS.SRTM1.(struc_name) = obsS_temp;
    end

    [file,path] = uigetfile('mod*ALTI3D*.shp','load file',[save_dir filesep 'modS_' char(demstr(3)) '.shp'],'MultiSelect','on');
    if isstr(file), file={file}; end
    for i=1:numel(file)
        struc_name = char(file(i)); struc_name = struc_name(1:end-4);
        modS_temp = shaperead([path char(file(i))]);
        modS.ALTI3D.(struc_name) = modS_temp;
        obsS_temp = shaperead([path strrep(char(file(i)),'modS','obsS')]);
        obsS.ALTI3D.(struc_name) = obsS_temp;
    end
    
    %plots
    field = 'length';
    climada_ls_probAssPlot(modS,obsS,field,'both',[0 2500],50,[0 2500])
    % climada_ls_probAssPlot(modS,obsS,'length','qq',[0 500])
    % climada_ls_probAssPlot(modS,obsS,'length','hist','',100,[0 2500])
    % climada_ls_probAssPlot(modS,obsS,'length','hist','',25,[0 1000])
end


%%
%histogram of distribution of phi for slides lenght = 0 and slides
%vs distribtuion of vmax phi for slides >0 --> maybe you see that it is
%just dependend on phi (or this is the case by definition to be clear)

%find zero_slides and plot phi dist values for zeros slides and non zero

zero_idx = find([modS_SRTM3.length]==0);
modS_zeros.SRTM3.SRTM3 = modS_SRTM3(zero_idx);
zero_idx = find([modS_SRTM1.length]==0);
modS_zeros.SRTM1.SRTM1 = modS_SRTM1(zero_idx);
zero_idx = find([modS_ALTI3D.length]==0);
modS_zeros.ALTI3D.ALTI3D = modS_ALTI3D(zero_idx);

%find slides bigger than zero
nonzero_idx = find([modS_SRTM3.length]>0);
modS_nonzeros.SRTM3.SRTM3 = modS_SRTM3(nonzero_idx);
nonzero_idx = find([modS_SRTM1.length]>0);
modS_nonzeros.SRTM1.SRTM1 = modS_SRTM1(nonzero_idx);
nonzero_idx = find([modS_ALTI3D.length]>0);
modS_nonzeros.ALTI3D.ALTI3D = modS_ALTI3D(nonzero_idx);

%dist zeros vs nonzeros phi
field = 'phi';
climada_ls_probAssPlot(modS_zeros,modS_nonzeros,field,'hist','',2.5,[0 60])
%dist zeros vs nonzeros vmax
field = 'vmax';
climada_ls_probAssPlot(modS_zeros,modS_nonzeros,field,'hist','',1,[1 11])

%%

%histogram of distribution of phi for slides lenght > threshold and slides
%vs distribtuion of vmax phi for slides >0 --> maybe you see that it is
%just dependend on phi (or this is the case by definition to be clear)

str_lgt_th = '1000'; %threshold for length
lgt_th = str2num(str_lgt_th);

%find slides smaller than threshold
idx = find([modS_SRTM3.length]<=lgt_th);
modS_stLgtTh.SRTM3.SRTM3 = modS_SRTM3(idx);
idx = find([modS_SRTM1.length]<=lgt_th);
modS_stLgtTh.SRTM1.SRTM1 = modS_SRTM1(idx);
idx = find([modS_ALTI3D.length]<=lgt_th);
modS_stLgtTh.ALTI3D.ALTI3D = modS_ALTI3D(idx);

%find slides bigger than threshold
idx = find([modS_SRTM3.length]>lgt_th);
modS_gtLgtTh.SRTM3.SRTM3 = modS_SRTM3(idx);
idx = find([modS_SRTM1.length]>lgt_th);
modS_gtLgtTh.SRTM1.SRTM1 = modS_SRTM1(idx);
idx = find([modS_ALTI3D.length]>lgt_th);
modS_gtLgtTh.ALTI3D.ALTI3D = modS_ALTI3D(idx);

%dist st threshold vs gt threshold for phi
field = 'phi';
climada_ls_probAssPlot(modS_stLgtTh,modS_gtLgtTh,field,'hist','',2.5,[0 60])
%find slides smaller than threshold
field = 'vmax';
climada_ls_probAssPlot(modS_stLgtTh,modS_gtLgtTh,field,'hist','',1,[1 11])

%%
%dist and qq when zero slides are not considered
%find slides bigger than zero
nonzero_idx = find([modS_SRTM3.length]>0);
modS_rmvzeros.SRTM3.SRTM3 = modS_SRTM3(nonzero_idx);
obsS_rmvzeros.SRTM3.SRTM3 = obsS_SRTM3(nonzero_idx);
nonzero_idx = find([modS_SRTM1.length]>0);
modS_rmvzeros.SRTM1.SRTM1 = modS_SRTM1(nonzero_idx);
obsS_rmvzeros.SRTM1.SRTM1 = obsS_SRTM1(nonzero_idx);
nonzero_idx = find([modS_ALTI3D.length]>0);
modS_rmvzeros.ALTI3D.ALTI3D = modS_ALTI3D(nonzero_idx);
obsS_rmvzeros.ALTI3D.ALTI3D = obsS_ALTI3D(nonzero_idx);

field = 'length';
climada_ls_probAssPlot(modS_rmvzeros,obsS_rmvzeros,'length','both',[0 2500],50,[0 2500])

%%
%get distance in dx and dy direction in observation --> translate dx and dy
%to dx dy it would have in smaller resolution --> get with pythagoras
%length it would have in smaller raster
deg_km = 111.32;
dlat_alti3d = abs(min(diff(lat_alti3d(:,1)))); 
dlon_alti3d = abs(min(diff(lon_alti3d(1,:))));
dy_alti3d = dlat_alti3d*(deg_km * 1000);
dx_alti3d = dlon_alti3d*cosd(mean(lat_alti3d(:,1)))*(deg_km * 1000);

dlat_srtm1 = abs(min(diff(lat_srtm1(:,1)))); 
dlon_srtm1 = abs(min(diff(lon_srtm1(1,:))));
dy_srtm1 = dlat_srtm1*(deg_km * 1000);
dx_srtm1 = dlon_srtm1*cosd(mean(lat_srtm1(:,1)))*(deg_km * 1000);

dlat_srtm3 = abs(min(diff(lat_srtm3(:,1)))); 
dlon_srtm3 = abs(min(diff(lon_srtm3(1,:))));
dy_srtm3 = dlat_srtm3*(deg_km * 1000);
dx_srtm3 = dlon_srtm3*cosd(mean(lat_srtm3(:,1)))*(deg_km * 1000);

obsS_norm.SRTM3.dL = [obsS_SRTM3.dL]; %not normalized horizontal length
obsS_norm.SRTM1.dL = [obsS_SRTM1.dL];
obsS_norm.ALTI3D.dL = [obsS_ALTI3D.dL];

%calculate observed length of model in x and y direction

%srtm3
obs_X = reshape([obsS_SRTM3.X],3,[]);
obs_Y = reshape([obsS_SRTM3.Y],3,[]);

x_diff = round(abs(obs_X(1,:)-obs_X(2,:))/dlon_srtm3)*dx_srtm3;
y_diff = round(abs(obs_Y(1,:)-obs_Y(2,:))/dlat_srtm3)*dy_srtm3;

dL = sqrt(x_diff.^2+y_diff.^2);

obsS_norm.SRTM3.SRTM3 = obsS_SRTM3;
obsS_norm.SRTM3.SRTM3_norm.dL = dL;

%srtm1
obs_X = reshape([obsS_SRTM1.X],3,[]);
obs_Y = reshape([obsS_SRTM1.Y],3,[]);

x_diff = round(abs(obs_X(1,:)-obs_X(2,:))/dlon_srtm1)*dx_srtm1;
y_diff = round(abs(obs_Y(1,:)-obs_Y(2,:))/dlat_srtm1)*dy_srtm1;

dL = sqrt(x_diff.^2+y_diff.^2);

obsS_norm.SRTM1.SRTM1 = obsS_SRTM1;
obsS_norm.SRTM1.SRTM1_norm.dL = dL;

%alti3d
obs_X = reshape([obsS_ALTI3D.X],3,[]);
obs_Y = reshape([obsS_ALTI3D.Y],3,[]);

x_diff = round(abs(obs_X(1,:)-obs_X(2,:))/dlon_alti3d)*dx_alti3d;
y_diff = round(abs(obs_Y(1,:)-obs_Y(2,:))/dlat_alti3d)*dy_alti3d;

dL = sqrt(x_diff.^2+y_diff.^2);

obsS_norm.ALTI3D.ALTI3D = obsS_ALTI3D;
obsS_norm.ALTI3D.ALTI3D_norm.dL = dL;

% %construct modS --> extract dL for plot
modS_norm.SRTM3.SRTM3_norm = modS_SRTM3;
modS_norm.SRTM1.SRTM1_norm = modS_SRTM1;
modS_norm.ALTI3D.ALTI3D_norm = modS_ALTI3D;
modS_norm.SRTM3.SRTM3 = modS_SRTM3;
modS_norm.SRTM1.SRTM1 = modS_SRTM1;
modS_norm.ALTI3D.ALTI3D = modS_ALTI3D;
% modS_norm.SRTM1.dL = [modS_SRTM1.dL];
% modS_norm.SRTM1.norm_dL = [modS_SRTM1.dL];
% modS_norm.ALTI3D.dL = [modS_ALTI3D.dL];
% modS_norm.ALTI3D.norm_dL = [modS_ALTI3D.dL];



climada_ls_probAssPlot(modS_norm,obsS_norm,'dL','qq',[0 2500],50,[0 2500])

%%
%plot cdf of with observed length and observed length when normalised
obsS_cdfNorm.SRTM3.SRTM3 = obsS_norm.SRTM3.SRTM3;
obsS_cdfNorm.SRTM3.SRTM3_Norm = obsS_norm.SRTM3.SRTM3_norm;
obsS_cdfNorm.SRTM1.SRTM1 = obsS_norm.SRTM1.SRTM1;
obsS_cdfNorm.SRTM1.SRTM1_Norm = obsS_norm.SRTM1.SRTM1_norm;
obsS_cdfNorm.ALTI3D.ALTI3D = obsS_norm.ALTI3D.ALTI3D;
obsS_cdfNorm.ALTI3D.ALTI3D_Norm = obsS_norm.ALTI3D.ALTI3D_norm;

field = 'dL';
climada_ls_probAssPlot(obsS_cdfNorm,obsS_cdfNorm,field,'both',[0 2500],50,[0 2500],1,[0 400])

%%

disp('hier')









end