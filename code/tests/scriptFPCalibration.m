function scriptFPCalibration(calculate,calculate2,calculate3)

% Script to calibrate flow path parameters
% INPUTS: 
%     calculate: 1 if area lenght start end of original high resolution Sample need to be calcualted
%     calculate2: 1 if snapping of starting area from high resolution to coarse resoution need to be caluculated
%     calculate3: 1 if area and length with different modelparameters need to be caculated
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180417, init

global climada_global
if ~climada_init_vars, return; end

%folder for calibration
dat_dir = 'C:\Users\Simon Rölli\Desktop\data\calibration';
field = 'ID';

if ~exist('calculate') calculate=0; end
if ~exist('calculate2') calculate2=0; end

%save directories
sav_dir_org = [dat_dir filesep 'orgData'];
sav_dir_snap = [dat_dir filesep 'snapData'];
sav_dir_cali = [dat_dir filesep 'caliData'];

%dem which shall be calibration high resolution (org) and low resolution (snap)
%orgSave = '2x3m';
file_dem_org = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat';
%calSave = '20x30m';
%file_dem_snap = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm1';%srtm1
file_dem_snap = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d';%alti3d ca. 7*10m

%%
%load data --> DEM org and snap/cali
%read in DEM (original alti3d with 2x3m resolution is best choice --> best results when assessing source
%area lenght and width of slides 
load(file_dem_org,'centroids');
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
orgLon = reshape(centroids.lon,n_lat,n_lon);
orgLat = reshape(centroids.lat,n_lat,n_lon);
orgElevation = reshape(centroids.elevation_m,n_lat,n_lon);

%get resolution of org DEM to define orgSave
deg_km = 111.32;
dlat = abs(min(diff(orgLat(:,1)))); 
dlon = abs(min(diff(orgLon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(orgLat(:,1)))*(deg_km * 1000); 
orgSave=[num2str(round(dx)) 'x' num2str(round(dy)) 'm'];

load(file_dem_snap,'centroids');
% calSave = '7x10m';
% load(,'centroids');
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

%get resolution of snap/cali DEM to define calSave
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
calSave=[num2str(floor(dx)) 'x' num2str(floor(dy)) 'm'];

clear centroids;
%%
%get from original polygon to original raster (with ca 2x3m resolution)
if calculate
    %read in shape file of past shallow landslides
    %S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');
    orgS = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');

    [orgS,orgPolyraster,orgStart,orgEnd] = climada_calibration_orgS2raster(orgLon,orgLat,orgElevation,orgS,field,0);

%     figure
%     plot([S.AREA_GIS],[S.P_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')
%     figure
%     plot([orgS.AREA_GIS],([orgS.R_AREA_noSL]-[orgS.AREA_GIS])./[orgS.AREA_GIS],'*')
%     xlabel('ArcGIS [m^2]')
%     ylabel('rel. Difference Matlab raster vs ArcGIS [%]')

    %calculate ratio of length/shape length
    C = num2cell([orgS.R_LGT_noSL]./[orgS.SHAPE_Leng]);
    [orgS.LGT_RAT_noSL] = C{:};
    %with slope
    C = num2cell([orgS.R_LGT_SL]./[orgS.SHAPE_Leng]);
    [orgS.LGT_RAT_SL] = C{:};

    shapewrite(orgS,[sav_dir_org filesep 'orgS_' orgSave '.shp']);
    grid2tiff(orgLon,orgLat,orgStart,[sav_dir_org filesep 'orgStart_' orgSave '.tif']);
    grid2tiff(orgLon,orgLat,orgEnd,[sav_dir_org filesep 'orgEnd_' orgSave '.tif']);
    grid2tiff(orgLon,orgLat,orgPolyraster,[sav_dir_org filesep 'orgPolyraster_' orgSave '.tif']);
end
%%
% snapping of high resolution (2x3m) starting points (org) in coarser resolution
% (snap)
if calculate2
   orgS = shaperead([sav_dir_org filesep 'orgS_' orgSave '.shp']);
   grid = GRIDobj([sav_dir_org filesep 'orgStart_' orgSave '.tif']);
   orgStart = flipud(grid.Z);
   grid = GRIDobj([sav_dir_org filesep 'orgEnd_' orgSave '.tif']);
   orgEnd = flipud(grid.Z);
   clear grid;
   
   %get length and starting/end points of slides when changing to coarser
   %resolution (snapping)
   [snapS,snapStart,snapEnd] = climada_calibration_snap(orgS,orgLon,orgLat,orgStart,orgEnd,lon,lat,elevation,field);
   
   shapewrite(snapS,[sav_dir_snap filesep 'snapS_' calSave '.shp']);
   grid2tiff(lon,lat,snapStart,[sav_dir_snap filesep 'snapStart_' calSave '.tif']);
   grid2tiff(lon,lat,snapEnd,[sav_dir_snap filesep 'snapEnd_' calSave '.tif']);
   
end


%%
%calculate slides parameter (area/length) for changing phi and vmax values,
%when trigger slides from source areas and model flow path
%calibrate flow path parameters (cali)

orgS = shaperead([sav_dir_org filesep 'orgS_' orgSave '.shp']);
% grid = GRIDobj(['C:\Users\Simon Rölli\Desktop\data\calibration\orgStart_' orgSave '.tif']);
% orgStart = flipud(grid.Z);
% grid = GRIDobj(['C:\Users\Simon Rölli\Desktop\data\calibration\orgEnd_' orgSave '.tif']);
% orgEnd = flipud(grid.Z);

snapS = shaperead([sav_dir_snap filesep 'snapS_' calSave '.shp']);
grid = GRIDobj([sav_dir_snap filesep 'snapStart_' calSave '.tif']);
snapStart = flipud(grid.Z);
% grid = GRIDobj([save_dir_snap filesep 'snapEnd_' calSave '.tif']);
% snapEnd = flipud(grid.Z);
clear grid;
  
%area of each raster cell --> considering slope
cell_area = climada_centroids_area(lon,lat,elevation,0);

%create folder to saveif not exising
sav_dir_cali = [sav_dir_cali filesep calSave];
if ~exist(sav_dir_cali,'dir'), mkdir(sav_dir_cali); end

if calculate3
    %%
    %create parameter sets
    sphi = 15;
    ephi = 30;
    svmax = 1;
    evmax = 12;

    tot_cal = (ephi-sphi+1)*(evmax-svmax+1);

    %constant parameters
    dH = repmat(0,1,tot_cal);
    exp = repmat(25,1,tot_cal);
    iT = repmat(0.0003,1,tot_cal);
    perWT = repmat([1 0.8 0.4 0 0 0 0.4 0.8],tot_cal,1);

    %to calibrate
    phi = sphi:1:ephi;
    vmax = svmax:1:evmax;
    [phi,vmax] = meshgrid(phi,vmax);
    phi = phi(:)';
    vmax = vmax(:)';
    %%
    %flow propagation with different parameters
    %iteration through calibration parameters; save S in the end
    for i=1:tot_cal
        %processing management
        fprintf('Parameter set %i of %i...\n\t',i,tot_cal);

        %extract flow path parameters
        flowPara.dH = dH(i);
        flowPara.exp = exp(i);
        flowPara.phi = phi(i);
        flowPara.vmax = vmax(i);
        flowPara.iT = iT(i);
        flowPara.perWT = perWT(i,:);

        %propagate and calculate length/area
        [S,~,~] = climada_ls_FPcalibration(lon,lat,elevation,flowPara,cell_area,snapStart,snapS,field);

        %remove point from iT and point and ' ' from perWT
        iT_str = num2str(iT(i));
        iT_str(iT_str=='.') = [];

        perWT_str = num2str(perWT(i,:));
        perWT_str(diff(double(perWT_str))==0) = [];
        perWT_str(perWT_str==' ') = '_';
        perWT_str(perWT_str=='.') = [];

        %create filename and save
        filename = [calSave '_phi' num2str(phi(i)) '_vmax' num2str(vmax(i)) '_exp' num2str(exp(i))...
            '_dH' num2str(dH(i)) '_iT' iT_str '_perWT' perWT_str];
        shapewrite(S,[sav_dir_cali filesep filename '.shp']);

    end

    fprintf(' done\n');
end


%%
%load propagated calibration datasets (with different parameters and write
%in a structure array
disp('')





end