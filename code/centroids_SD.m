function centroids=centroids_SD(centroids, check_plots)
% Assign soil depth (SD) values to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_SD
% PURPOSE:
%   Determine soil depth (SD) for given centroids by
%   reading a grayscale GeoTIFF of global annual mean SD.
%   The global average annual SD GeoTIFF image derived from http://soilgrids.org/
%   
%   It needs to be placed in data/system of the climada module flood
%   The code automatically reads the GeoTIFF from that URL if it does not
%   exist in the data folder of the flood module (but this is VERY slow...).
%   soil depth values are in m
% CALLING SEQUENCE:
%   centroids = centroids_SD(centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_SD(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMSDERS:
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, SD is calculated even if the centroids
%   already have a field 'soil depth' (default is 0)
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'SD_mm' containing info about soil depth to bedrock in mm
% NOTE:
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmaail.com, 20150414
%-

% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids',  'var') || isempty(centroids),   climada_centroids_load; end
if ~exist('check_plots','var') || isempty(check_plots), check_plots = 0;        end

% locate the module's data folder
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% PARAMETERS
%
% the file with the SD data
SD_file = [module_data_dir filesep 'system' filesep 'soil_data' filesep 'BDRICM_02_apr_2014.tif'];
%
% URL where the SD file can be downloaded
SD_file_URL  = 'ftp://soilgrids:soilgrids@ftp.soilgrids.org/data/recent/BDRICM_02_apr_2014.tif.gz';
%
%
% north south extension of the SD image (longitude is covered in the full
% range from -180 to 180)
min_lat=-90; % degree
max_lat= 90;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global SD map)
bbox = [min(centroids.lon), min(centroids.lat),max(centroids.lon), max(centroids.lat)];

[fP,fN,~] = fileparts(SD_file);
SD_file_mat = [fP filesep fN '.mat'];
if exist(SD_file_mat,'file') %climada_check_matfile(SD_file,SD_file_mat)
    % .mat file exists, we just have to load it
    load(SD_file_mat);
elseif exist(SD_file,'file')
    % GeoTIFF file exists; read it
    fprintf('reading soil depth data from %s\n',SD_file)
    full_SD_img = imread(SD_file);
    fprintf('saving as .mat file...')
    save(SD_file_mat,'full_SD_img');
    fprintf('done \n')
else
    % Read the GeoTIFF file from the web
    fprintf('reading soil depth data from web:%s\n',SD_file_URL)
    gunzip(SD_file_URL,fP);
    full_SD_img = imread(SD_file);
    fprintf('saving as .mat file...')
    save(SD_file_mat,'full_SD_img');
    fprintf('done \n')
end

full_SD_img=flipud(full_SD_img); % switch for correct order in latitude (images are saved 'upside down')

% Crop the image to the centroids' bounding box (for speedup)
lon         = 360*(1:size(full_SD_img,2))/size(full_SD_img,2)+(-180); % -180..180
lat         = (max_lat-min_lat)*(1:size(full_SD_img,1))/size(full_SD_img,1)+min_lat;
lon_crop_ndx= lon>=bbox(1) & lon<=bbox(3);
lat_crop_ndx= lat>=bbox(2) & lat<=bbox(4);
img         = full_SD_img(lat_crop_ndx,lon_crop_ndx);
lon         = lon(lon_crop_ndx);
lat         = lat(lat_crop_ndx);
[LON,LAT]	= meshgrid(lon,lat); % construct regular grid

% depth to bedrock up to a max of 240cm (see original documentation), hence
% set all vals >240 to 0 and convert to mm
img = double(img);
img(img>240) = 0;
img = img .* 10;

fprintf('assigning soil depth values to centroids... ')
centroids.SD_mm= interp2(LON,LAT,img,centroids.lon,centroids.lat,'cubic');
fprintf(' done\n');

if check_plots
    % convert to double (from uint16)
    soil_depth=griddata(centroids.lon,centroids.lat,centroids.SD_mm,LON,LAT);%  double(img);
    % plot the image (kind of 'georeferenced')
    figure('color', 'w')
    hold on
    climada_plot_world_borders
    axis([bbox(1) bbox(3) bbox(2) bbox(4)])
    pcolor(LON,LAT,soil_depth); 
    colorbar
    shading flat
    if isfield(centroids,'country_name')
        title_string = sprintf('soil depth (mm), %s', centroids.country_name{1});
    else
        title_string = 'soil depth (mm/day)';
    end
    title(title_string)
    hold off
end