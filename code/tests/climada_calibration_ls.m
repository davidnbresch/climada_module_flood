function climada_calibration_ls(org_centroids,cal_centroids,orgS,flowPara,sav_dir,subS,snapS)

% Script to calibrate flow path parameters 
% MODULE:
%   flood
% NAME:
%   climada_calibration_ls
% PURPOSE:
%   Is used to calibrate 
% CALLING SEQUENCE:
%   climada_calibration_ls(calculate,calculate2,calculate3)
% EXAMPLE:
%   
% INPUTS:
%   org_centroids:  original/high resolution DEM (as centroid)
%                   transformation from polygon to raster will be based on
%                   this resolution. The area/length and other scores for
%                   the calibration will be derived at this DEM-resolution
%   cal_centroids:  snapping and calibration/low resolution DEM (as
%                   centroid). Snapping is based on this resolution/DEM.
%                   Scores are caculated from this resolution and starting
%                   points (nearest grid points to starting point in
%                   original DEM) are latter used as sources cells for
%                   propagation model.
%   orgS:           original shape file with polygons representing shallow
%                   landslides. Need to include .X (vector with lon 
%                   coordinates of  vertices of polygons), .Y (vector 
%                   with lat coordinates vertices of polygons) and .ID (numerical
%                   value to label slides) fields
%   save_dir:       string with path which defines directory to save files
%                   If not already existing, program will greate folder for
%                   original, snapping and calibration data in the
%                   save-directory.
% OPTIONAL INPUT PARAMETERS:
%   snapS:      Provides structure of slides after    
%   calculate:  1 (default=0) if area lenght start end of original high 
%               resolution Sample need to be calcualted (When transforming
%               polygones to a raster file).
%   calculate2: 1 (default=0) if snapping of starting area from high 
%               resolution to coarse resoution need to be caluculated
%   calculate3: 1 (default=0) if area and length with different
%               model parameters need to be caculated
% OUTPUTS:
%   no output but saves shapefile of original raster/snapped and calibrated
%   /modelled shape files in defined savepaths
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180417, init
% Thomas Rölli, thomasroelli@gmail.com, 20180522, rename and restructure

global climada_global
if ~climada_init_vars, return; end

%folder for calibration
%dat_dir = 'C:\Users\Simon Rölli\Desktop\data\calibration';
field = 'ID';

if ~exist('org_centroids') return; end
if ~exist('cal_centroids') return; end
if ~exist('orgS') return; end
if ~exist('flowPara') return; end
if ~exist('sav_dir') return; end
if ~exist('subS','var') subS=[]; end
if ~exist('snapS') snapS=[]; end
    
%create directory to save if not existing
if ~exist(sav_dir,'dir'), mkdir(sav_dir); end 

%folder for calibration
%dat_dir = 'C:\Users\Simon Rölli\Desktop\data\calibration';
%slide labelling which is used and assinged in the follow shape files
field = 'ID';

%save directories
sav_dir_org = [sav_dir filesep 'orgData'];
sav_dir_sub = [sav_dir filesep 'subData'];
sav_dir_snap = [sav_dir filesep 'snapData'];
sav_dir_cali = [sav_dir filesep 'caliData'];
if ~exist(sav_dir_sub,'dir'), mkdir(sav_dir_sub); end
if ~exist(sav_dir_snap,'dir'), mkdir(sav_dir_snap); end
if ~exist(sav_dir_cali,'dir'), mkdir(sav_dir_cali); end

%%
%read in DEM (original alti3d with 2x3m resolution is best choice --> best results when assessing source
%area lenght and width of slides 
centroids = org_centroids;
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

%DEM snap/cali
centroids = cal_centroids;
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
if isempty(subS)
    %read in shape file of past shallow landslides
    %S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');
    %orgS = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');

    [subS,~,~,~] = climada_calibration_orgS2raster(orgLon,orgLat,orgElevation,orgS,field,0);

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
%     C = num2cell([orgS.R_LGT_noSL]./[orgS.SHAPE_Leng]);
%     [orgS.LGT_RAT_noSL] = C{:};
%     %with slope
%     C = num2cell([orgS.R_LGT_SL]./[orgS.SHAPE_Leng]);
%     [orgS.LGT_RAT_SL] = C{:};

    shapewrite(orgS,[sav_dir_org filesep 'orgS_' orgSave '.shp']);
    shapewrite(subS,[sav_dir_sub filesep 'subS_' orgSave '.shp']);
%     grid2tiff(orgLon,orgLat,orgStart,[sav_dir_org filesep 'orgStart_' orgSave '.tif']);
%     grid2tiff(orgLon,orgLat,orgEnd,[sav_dir_org filesep 'orgEnd_' orgSave '.tif']);
%     grid2tiff(orgLon,orgLat,orgPolyraster,[sav_dir_org filesep 'orgPolyraster_' orgSave '.tif']);
end
%%
% snapping of high resolution (2x3m) starting points (org) in coarser resolution
% (snap)
if isempty(snapS)
   %subS = shaperead([sav_dir_sub filesep 'subS_' orgSave '.shp']);
%    grid = GRIDobj([sav_dir_org filesep 'orgStart_' orgSave '.tif']);
%    orgStart = flipud(grid.Z);
%    grid = GRIDobj([sav_dir_org filesep 'orgEnd_' orgSave '.tif']);
%    orgEnd = flipud(grid.Z);
%    clear grid;
   
   %get length and starting/end points of slides when changing to coarser
   %resolution (snapping)
   [snapS,~,~] = climada_calibration_snap(subS,orgLon,orgLat,lon,lat,elevation,field);
   
   shapewrite(snapS,[sav_dir_snap filesep 'snapS_' calSave '.shp']);
%    grid2tiff(lon,lat,snapStart,[sav_dir_snap filesep 'snapStart_' calSave '.tif']);
%    grid2tiff(lon,lat,snapEnd,[sav_dir_snap filesep 'snapEnd_' calSave '.tif']);
   
end


%%
%calculate slides parameter (area/length) for changing phi and vmax values,
%when trigger slides from source areas and model flow path
%calibrate flow path parameters (cali)
  
%area of each raster cell --> considering slope
cell_area = climada_centroids_area(lon,lat,elevation,0);

%create folder to saveif not exising
sav_dir_cali = [sav_dir_cali filesep calSave];
if ~exist(sav_dir_cali,'dir'), mkdir(sav_dir_cali); end

   
%flow propagation with different parameters
%iteration through calibration parameters; save S in the end
for i=1:numel(flowPara.phi)
    %processing management
    fprintf('Parameter set %i of %i...\n\t',i,numel(flowPara.phi));

    %extract flow path parameters
    subFlowPara.dH = flowPara.dH(i);
    subFlowPara.exp = flowPara.exp(i);
    subFlowPara.phi = flowPara.phi(i);
    subFlowPara.vmax = flowPara.vmax(i);
    subFlowPara.iT = flowPara.iT(i);
    subFlowPara.perWT = flowPara.perWT(i,:);

    %propagate and calculate length/area
    [S,~,~] = climada_calibration_propagate(lon,lat,elevation,subFlowPara,cell_area,snapS,field);

    %remove point from iT and point and ' ' from perWT
    iT_str = num2str(subFlowPara.iT);
    iT_str(iT_str=='.') = [];

    perWT_str = num2str(subFlowPara.perWT);
    perWT_str(diff(double(perWT_str))==0) = [];
    perWT_str(perWT_str==' ') = '_';
    perWT_str(perWT_str=='.') = [];

    %create filename and save
    filename = [calSave '_phi' num2str(subFlowPara.phi) '_vmax'...
        num2str(subFlowPara.vmax) '_exp' num2str(subFlowPara.exp)...
        '_dH' num2str(subFlowPara.dH) '_iT' iT_str '_perWT' perWT_str];
    shapewrite(S,[sav_dir_cali filesep filename '.shp']);
    
end

fprintf(' done\n');

    


%%
%load propagated calibration datasets (with different parameters and write
%in a structure array

% files = dir([sav_dir_cali '/*.shp']);
% res_cali = repmat(struct('resol',[]),1,numel(files));
% 
% for i=1:numel(files)
%    fileparts = strsplit(files(i).name,'_');
%    res_cali(i).resol = fileparts(1);
%    %phi
%    num_idx = regexp(fileparts{2},'\d');
%    res_cali(i).phi = str2double(extractAfter(fileparts{2},num_idx(1)-1));
%    %vmax
%    num_idx = regexp(fileparts{3},'\d');
%    res_cali(i).vmax = str2double(extractAfter(fileparts{3},num_idx(1)-1));
%    %exp
%    num_idx = regexp(fileparts{4},'\d');
%    res_cali(i).exp = str2double(extractAfter(fileparts{4},num_idx(1)-1));
%    %dH
%    num_idx = regexp(fileparts{5},'\d');
%    res_cali(i).dH = str2double(extractAfter(fileparts{5},num_idx(1)-1));
%    %iT
%    num_idx = regexp(fileparts{6},'\d');
%    tmp_str = extractAfter(fileparts{6},num_idx(1)-1);
%    res_cali(i).iT = str2double([tmp_str(1) '.' tmp_str(2:end)]);
%    %perWTS
%    num_idx = regexp([fileparts{7:end}],'\d');
%    tmp_str = extractBetween([fileparts{7:end}],num_idx(1),num_idx(end));
%    res_cali(i).perWT = tmp_str;
%    
%    %open corresponding caliS
%    caliS = shaperead([files(i).folder filesep files(i).name]);
%    
%    obs = [snapS.length];%observed lengths
%    pred = [caliS.length];%predicted lengths
%    
%    res_cali(i).rmse = sqrt(mean((obs-pred).^2,'omitnan'));
%    
%    res_cali(i).source = [files(i).folder filesep files(i).name];
% end
% 
% %plot RMSE vs phi for different vmax
% if 0
%     vmax = [res_cali.vmax];
%     vmax_unique = unique(vmax);
%     phi = [res_cali.phi];
%     phi_unique = unique(phi);
%     rmse = [res_cali.rmse];
%     figure
%     for i=1:numel(vmax_unique)
%         idx_vmax = find(vmax == vmax_unique(i));
%         subplot(3,4,i)
%         plot(phi(idx_vmax),rmse(idx_vmax),'-*')
%         title(['vmax = ' num2str(vmax_unique(i))])
%         xlim([min(phi_unique) max(phi_unique)])
%         ylim([0 500])
%         xlabel('phi')
%         ylabel('RMSE')
%     end
% 
% %     figure
% %     for i=1:numel(phi_unique)
% %         idx_phi = find(phi == phi_unique(i));
% %         subplot(4,4,i)
% %         plot(vmax(idx_phi),rmse(idx_phi),'-*')
% %         title(['phi = ' num2str(phi_unique(i))])
% %         xlim([min(vmax_unique) max(vmax_unique)])
% %         ylim([0 500])
% %         xlabel('vmax')
% %         ylabel('RMSE')
% %     end
% end
% 
% %%
% %removement of some slides 
% %%%%%%%%%%%%%%%
% %%%%%ToDo%%%%%%
% %write function where thresholds can be set again
% %plot modelled vs original with and without removed slides
% %calculate difference (evt. RMSE) of each slide and look at characteristics
% % of them. e.g. points which are far away from 1:1 line
% %evt boxplots of lenght for each parameter set
% %maybe set for slide without propagation not lenght = 0 but with e.g.
% %max(dy,dx)/2
% %
% %done: 
% %remove too long and too large slides --> 
% %remove too small slides (area and/or length) according to unitarea of
% % resolution --> done
% %remove slide at steep catchment but which are too small --> we want worst
% %case --> fully saturated, model not able to resolve short slide at steep
% % catchments (therefore plot slop at source vs length)
% %remove slide at flat areas (low slope) --> model not able to model them 
% %maybe assign length value different from 0 for slides which do not
% % propagate
% %plots of lines which have color according to e.g. slope (maybe in ArcGIS)
% %Problems: at moment best parameters --> highest phi and lowest vmax --> a
% % lot of slides do not propagate because local maximum slope gradient is
% % smaller than phi --> lenght = 0 for all this slides --> we get better
% % RMSE --> maybe remove some of this small values
% %
% 
% %save slides which shall be removed from further considerations (1) and set
% %thresholds
% rmv = zeros(size([orgSubS.(field)])); %all removed slides
% 
% rmv_S2R = zeros(size([orgSubS.(field)])); %already removed in S2R or snap
% 
% rmv_large = zeros(size([orgSubS.(field)])); %too large (area)
% max_area_th = 30000; %[m^2] according to Millage etal area (10^1–10^4 m2).
% max_lgt_th = 1000;
% 
% rmv_small_area = zeros(size([orgSubS.(field)])); %too small (area) compared to unitarea
% min_area_th = dx*dy;
% rmv_small_lgt = zeros(size([orgSubS.(field)])); %too small (length) 
% min_lgt_th = 100;
% 
% rmv_minslope = zeros(size([orgSubS.(field)])); %too small slopes at start --> will not propagate
% min_slope_th = 15; %degrees minimum
% 
% rmv_slVSlgt = zeros(size([orgSubS.(field)])); %too small length for steep slope --> model will be too long
% slVSlgt_th = [25 80]; %slides with start slope over th(1) and lenght under (2)
% 
% %extract data
% area = [orgSubS.area];
% length = [orgSubS.length];
% st_slope = [orgSubS.st_slope];
% sl_slope = [orgSubS.sl_slope]; %slope of slide (line from start to end)
% max_st_slope = [snapS.max_srcslop]; %slope gradient maximum considert 8 neighbours
% 
% %%%%%%%%%%%%%%%
% %starting plots (nothing removed)
% % figure
% % plot(area,st_slope,'*'); xlabel('area'); ylabel('st_slope')
% % figure
% % plot(area,sl_slope,'*'); xlabel('area'); ylabel('sl_slope')
% % figure
% % plot(length,sl_slope,'*'); xlabel('length'); ylabel('sl_slope')
% % figure
% % plot(length,st_slope,'*'); xlabel('length'); ylabel('st_slope') %we dont see any correlation
% %%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%
% %plot([orgSubS.length],[snapS.length],'*')
% 
% % plot([orgSubS.area],[orgS.R_AREA_SL],'*') %--> 1:1 line
% % plot([orgSubS.length],[orgS.R_LGT_SL],'*') % --> 1:1 line
% 
% %%%%%%%%%%%%
% %save slide which have already been removed in orgS (2) or snapS (1)
% rmv_S2R([snapS.removed]~=0) = 1;
% %%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%
% %remove slides which are too long/big
% %plot([orgSubS.length]) %lenght threshold around 500m
% %plot([orgSubS.area]) %area threshold around 30000m^2
% rmv_large(area>max_area_th) = 1;
% rmv_large(length>max_lgt_th) = 1;
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %remove slides which are too small
% rmv_small_area(area<(min_area_th)) = 1;
% 
% rmv_small_lgt(length<(min_lgt_th)) = 1;
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %remove slides which have a too small start slope
% rmv_minslope(st_slope<min_slope_th) = 1;
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %remove slides which are too short for steep slope
% rmv_slVSlgt(length<slVSlgt_th(2) & st_slope>slVSlgt_th(1)) = 1;
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %set together rmv
% rmv = rmv_S2R+rmv_large+rmv_small_area+rmv_small_lgt+rmv_minslope+rmv_slVSlgt;
% rmv = logical(rmv);
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %caculate RMSE again without slides in rmv
% res_cali_rmv = repmat(struct('resol',[]),1,numel(files));
% for i=1:numel(files)
%    fileparts = strsplit(files(i).name,'_');
%    res_cali_rmv(i).resol = fileparts(1);
%    %phi
%    num_idx = regexp(fileparts{2},'\d');
%    res_cali_rmv(i).phi = str2double(extractAfter(fileparts{2},num_idx(1)-1));
%    %vmax
%    num_idx = regexp(fileparts{3},'\d');
%    res_cali_rmv(i).vmax = str2double(extractAfter(fileparts{3},num_idx(1)-1));
%    %exp
%    num_idx = regexp(fileparts{4},'\d');
%    res_cali_rmv(i).exp = str2double(extractAfter(fileparts{4},num_idx(1)-1));
%    %dH
%    num_idx = regexp(fileparts{5},'\d');
%    res_cali_rmv(i).dH = str2double(extractAfter(fileparts{5},num_idx(1)-1));
%    %iT
%    num_idx = regexp(fileparts{6},'\d');
%    tmp_str = extractAfter(fileparts{6},num_idx(1)-1);
%    res_cali_rmv(i).iT = str2double([tmp_str(1) '.' tmp_str(2:end)]);
%    %perWTS
%    num_idx = regexp([fileparts{7:end}],'\d');
%    tmp_str = extractBetween([fileparts{7:end}],num_idx(1),num_idx(end));
%    res_cali_rmv(i).perWT = tmp_str;
%    
%    %open corresponding caliS
%    caliS = shaperead([files(i).folder filesep files(i).name]);
%    
%    obs = [snapS.length];%observed lengths
%    pred = [caliS.length];%predicted lengths
%    
%    %set rmv to nan
%    obs(rmv) = nan;
%    pred(rmv) = nan;
%    
%    %calculate RMSE
%    res_cali_rmv(i).rmse = sqrt(mean((obs-pred).^2,'omitnan'));
%    
%    %save source file
%    res_cali_rmv(i).source = [files(i).folder filesep files(i).name];
% end
% %%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%
% %plot RMSE vs phi for rmv
% if 0
%     vmax = [res_cali_rmv.vmax];
%     vmax_unique = unique(vmax);
%     phi = [res_cali_rmv.phi];
%     phi_unique = unique(phi);
%     rmse = [res_cali_rmv.rmse];
%     figure
%     for i=1:numel(vmax_unique)
%         idx_vmax = find(vmax == vmax_unique(i));
%         subplot(3,4,i)
%         plot(phi(idx_vmax),rmse(idx_vmax),'-*')
%         title(['vmax = ' num2str(vmax_unique(i))])
%         xlim([min(phi_unique) max(phi_unique)])
%         ylim([0 500])
%         xlabel('phi')
%         ylabel('RMSE')
%     end
% 
% %     figure
% %     for i=1:numel(phi_unique)
% %         idx_phi = find(phi == phi_unique(i));
% %         subplot(4,4,i)
% %         plot(vmax(idx_phi),rmse(idx_phi),'-*')
% %         title(['phi = ' num2str(phi_unique(i))])
% %         xlim([min(vmax_unique) max(vmax_unique)])
% %         ylim([0 500])
% %         xlabel('vmax')
% %         ylabel('RMSE')
% %     end
% end
% 
% %%
% %%%%%%%%%%%%
% %make comparisons modelled vs original
% %
% %read in data
% caliS = shaperead(res_cali(151).source); %vmax=4, phi=27 we have local minimum in RMSE
% mod_length = [caliS.length];
% mod_area = [caliS.area];
% 
% %plot slides at slope below 
% min_slope_th = 27; %degrees minimum
% figure
% subplot(1,2,1)
% plot(length,mod_length,'*')
% title(['slides with start max gradient<' num2str(min_slope_th) ' [degrees]'])
% xlabel('observed [m]'); ylabel('modelled [m]'); xlim([0 1000]); ylim([0 1000]);
% hold on
% plot(length(max_st_slope<min_slope_th),mod_length(max_st_slope<min_slope_th),'*','color','red')
% subplot(1,2,2)
% plot(area,mod_area,'*')
% xlabel('observed [m^2]'); ylabel('modelled [m^2]'); xlim([0 30000]); ylim([0 30000]);
% hold on
% plot(area(st_slope<min_slope_th),mod_area(st_slope<min_slope_th),'*','color','red')
% 
% %other things to plot
% max_area_th = 30000; %[m^2] according to Millage etal area (10^1–10^4 m2).
% max_lgt_th = 1000;
% min_area_th = dx*dy;
% min_lgt_th = 100;
% slVSlgt_th = [25 80]; %slides with start slope over th(1) and length under (2)
% 
% %%%%%%%%%
% %plot relative difference vs slope....
% 
% reldiff_lgt = (mod_length-length)./length



end