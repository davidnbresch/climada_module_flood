function centroids=centroids_WHC(centroids, check_plots)
% Assign available water holding capacity (WHC) values to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_WHC
% PURPOSE:
%   Determine water holding capacity (WHC) for given centroids by
%   reading a grayscale GeoTIFF of global annual mean WHC.
%   The available water-holding capacity describes the amount of water
%   available for evapotranspiration after drainage and is defined as the
%   field capacity minus the permanent wilting point.
%   The global mean annual FC GeoTIFF image is part of the Global Gridded
%   Surfaces of Selected Soil Characteristics data set, which was developed
%   by the Global Soil Data Task of the International Geosphere-Biosphere
%   Programme. It can be downloaded here:
%   http://webmap.ornl.gov/wcsdown/wcsdown.jsp?dg_id=569_3
%   The image needs to be placed in data/system of the climada module
%   flood.
%   WHC values are in mm.
%   For more information on the Global Soil Data Task, see
%   https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
% CALLING SEQUENCE:
%   centroids = centroids_WHC(centroids, check_plots, force_recalc)
% EXAMPLE:
%   centroids = centroids_WHC(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, FC is calculated even if the centroids
%   already have a field 'water_holding_capacity' (default is 0)
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'WHC_mm' (water holding capacity), which contains the annual mean available
%   water holding capacity
%   (in mm) of each centroid
% NOTE:
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150325, initial
% Gilles Stassen, gillesstassen@hotmail.com, 20150423, cleanup
%-


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids), climada_centroids_load; end
if ~exist('check_plots', 'var') || isempty(check_plots), check_plots = 0;   end

% locate the module's data folder
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% PARAMETERS
%
% the file with the WHC data
WHC_file=[module_data_dir filesep 'system' filesep 'wcs_WHC.tiff'];
%
% Missing data value, to be found in the documentation here:
% ftp://www.daac.ornl.gov/data/global_soil/IGBP-SurfaceProducts/comp/readme.txt
missing_data_value = -2;
%
% north south range of the FC image (longitude is covered in the full
% range from -180 to 180)
min_lat=-56.5; % degree
max_lat= 90;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global ET map)
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon), min(centroids.lat), max(centroids.lon), max(centroids.lat)];

[fP,fN,~] = fileparts(WHC_file);
WHC_file_mat = [fP filesep fN '.mat'];
if climada_check_matfile(WHC_file,WHC_file_mat)
    % .mat file exists, we just have to load it
    load(WHC_file_mat);
elseif exist(WHC_file,'file')
    % GeoTIFF file exists; read it
    fprintf('Reading water holding capacity data from %s\n',WHC_file)
    full_WHC_img = imread(WHC_file);
    save(WHC_file_mat,'full_WHC_img');
else
    % Error message
    cprintf([206 50 50]/255,['ERROR: missing tiff image %s - '...
        'download data here: http://webmap.ornl.gov/wcsdown/wcsdown.jsp?dg_id=569_3'],WHC_file);
    return
end

full_WHC_img=full_WHC_img(end:-1:1,:); % switch for correct order in latitude (images are saved 'upside down')

% Crop the image to the centroids' bounding box (for speedup)
lon         = 360*(1:size(full_WHC_img,2))/size(full_WHC_img,2)+(-180); % -180..180
lat         = (max_lat-min_lat)*(1:size(full_WHC_img,1))/size(full_WHC_img,1)+min_lat;

bb          = max(max(diff(lon)),max(diff(lat))); %buffer

lon_crop_ndx= lon >= bbox(1)-bb & lon <= bbox(3)+bb;
lat_crop_ndx= lat >= bbox(2)-bb & lat <= bbox(4)+bb;
img         = full_WHC_img(lat_crop_ndx,lon_crop_ndx);
lon         = lon(lon_crop_ndx);
lat         = lat(lat_crop_ndx);
[LON,LAT]   = meshgrid(lon,lat); % construct regular grid

% % Set zeros to a fill value (defined by the median of all positive
% % data points in the section of the map we are looking at)
% img(img==0) = median(median(img(img>0)));
% % Set the NaNs (i.e. waterborne pixels) to zero
% img(img==missing_data_value) = 0;

fprintf('assigning water holding capacity values to centroids... ')
centroids.WHC_mm= interp2(LON,LAT,double(img),centroids.lon,centroids.lat,'linear');
fprintf(' done \n');

if check_plots
    % convert to double
    water_holding_cap=double(img);
    % plot the image (kind of 'georeferenced')
    pcolor(LON,LAT,water_holding_cap); colorbar
    shading flat
    if isfield(centroids,'admin0_name')
        title_string = sprintf(['Available water holding capacity'...
            '(mm), %s'], centroids.admin0_name);
    else
        title_string = 'Available water holding capacity (mm)';
    end
    title(title_string)
    hold on
    set(gcf,'Color',[1 1 1]) % white figure background
end