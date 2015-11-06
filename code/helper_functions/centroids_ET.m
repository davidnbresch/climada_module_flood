function centroids=centroids_ET(centroids, check_plots)
% Assign Evapotranspiration (ET) values to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_ET
% PURPOSE:
%   Determine Evapotranspiration (ET) for given centroids by
%   reading a grayscale GeoTIFF of global annual mean ET.
%   The global average annual ET GeoTIFF image derived from 2010-2013 MODIS
%   data can be downloaded here:
%   ftp://ftp.ntsg.umt.edu/pub/MODIS/NTSG_Products/MOD16/MOD16A3.105_MERRAGMAO/Geotiff/MOD16A3_ET_2000_to_2013_mean.tif
%   It needs to be placed in data/system of the climada module flood
%   The code automatically reads the GeoTIFF from that URL if it does not
%   exist in the data folder of the flood module (but this is VERY slow...).
%   Evapotranspiration values are in mm/day.
% CALLING SEQUENCE:
%   centroids = centroids_ET(centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_ET(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, ET is calculated even if the centroids
%   already have a field 'evapotranspiration' (default is 0)
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'ET_mm_day', which contains the annual mean evapotranspiration
%   (in mm/day) for each centroid
% NOTE:
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150320, initial
% Gilles Stassen, gillesstassen@hotmail.com, 20150407, cleanup, units mm/yr -> mm/day
% Gilles Stassen, gillesstassen@hotmaail.com, 20150414, speedup
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
% the file with the ET data
ET_file=[module_data_dir filesep 'system' filesep 'MOD16A3_ET_2000_to_2013_mean.tif'];
%
% URL where the ET file can be downloaded
ET_file_URL = ['ftp://ftp.ntsg.umt.edu/pub/MODIS/NTSG_Products/MOD16/MOD16A3.105_MERRAGMAO/Geotiff/MOD16A3_ET_2000_to_2013_mean.tif'];
%
% Some values have a special meaning; see documentation of the MODIS Global
% Evapotranspiration Project (MOD16) on
% http://www.ntsg.umt.edu/project/mod16
% Note: We will set all of them to 0, but with the information below, one
% could also treat different 'special cases' separately
% For the annual ET, the valid value range is 0-65500.
fill_value = 65535;     % Fill value, out of the earth
water_body = 65534;     % Water body
barren = 65533;         % Barren or sparsely vegetated
snow = 65532;           % Permanent snow and ice
wetland = 65531;        % Permanent wetland
urban = 65530;          % Urban or Built-up
unclassified = 65529;   % Unclassified
%
% north south extension of the ET image (longitude is covered in the full
% range from -180 to 180)
min_lat=-60; % degree
max_lat= 80;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global ET map)
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon)-1, min(centroids.lat)-1,max(centroids.lon)+1, max(centroids.lat)+1];

[fP,fN,~] = fileparts(ET_file);
ET_file_mat = [fP filesep fN '.mat'];
if exist(ET_file_mat,'file') %climada_check_matfile(ET_file,ET_file_mat)
    % .mat file exists, we just have to load it
    load(ET_file_mat);
elseif exist(ET_file,'file')
    % GeoTIFF file exists; read it
    fprintf('reading evapotranspiration data from %s\n',ET_file)
    full_ET_img = imread(ET_file);
    fprintf('saving as .mat file...')
    save(ET_file_mat,'full_ET_img');
    fprintf('done \n')
else
    % Read the GeoTIFF file from the web
    fprintf('reading evapotranspiration data from web:%s\n',ET_file_URL)
    full_ET_img = imread(ET_file_URL);
    fprintf('saving as .mat file...')
    save(ET_file_mat,'full_ET_img');
    fprintf('done \n')
end

full_ET_img=full_ET_img(end:-1:1,:); % switch for correct order in latitude (images are saved 'upside down')

% Set all 'special cases' defined in the parameters section to 0 (i.e.,
% to black)
% full_ET_img(full_ET_img>65500) = 0;

% Crop the image to the centroids' bounding box (for speedup)
lon         = 360*(1:size(full_ET_img,2))/size(full_ET_img,2)+(-180); % -180..180
lat         = (max_lat-min_lat)*(1:size(full_ET_img,1))/size(full_ET_img,1)+min_lat;
lon_crop_ndx= lon>=bbox(1) & lon<=bbox(3);
lat_crop_ndx= lat>=bbox(2) & lat<=bbox(4);
img         = full_ET_img(lat_crop_ndx,lon_crop_ndx); img(img>65500) = nanmean(nanmean(img))./2;
lon         = lon(lon_crop_ndx);
lat         = lat(lat_crop_ndx);
[LON,LAT]	= meshgrid(lon,lat); % construct regular grid

fprintf('assigning evapotranspiration values to centroids... ')
centroids.ET_mm_day= interp2(LON,LAT,double(img),centroids.lon,centroids.lat,'linear');


% perhaps a better method fro .onLand
% centroids.onLand(centroids.ET_mm_day == 65534)      = 0;
centroids.ET_mm_day(centroids.ET_mm_day > 65500)    = NaN;
centroids.ET_mm_day = centroids.ET_mm_day./65500;
% % Determine the pixel closest to each centroid
% n_centroids = length(centroids.centroid_ID);
% X_1D = reshape(LON,[1,numel(LON)]);
% Y_1D = reshape(LAT,[1,numel(LAT)]);
% fprintf('assigning evapotranspiration values to centroids... ')
% for centroid_i=1:n_centroids
%     distances=climada_geo_distance(centroids.lon(centroid_i),centroids.lat(centroid_i),X_1D,Y_1D);
%     [~,min_dist_index] = min(distances);
%     [img_row, img_col] = ind2sub(size(LON),min_dist_index);
%     centroids.ET_mm_day(centroid_i) = double(img(img_row, img_col));
% end

% convert range of greyscale values (uint16, i.e. values range from 0
% to 65535) into actual evapotranspiration (in mm/day). (the factor 0.1
% is taken from MODIS website)
centroids.ET_mm_day = centroids.ET_mm_day./0.1;
fprintf(' done\n');

if check_plots
    % convert to double (from uint16)
    evapo=griddata(centroids.lon,centroids.lat,centroids.ET_mm_day,LON,LAT);%  double(img);
    % plot the image (kind of 'georeferenced')
    pcolor(LON,LAT,evapo); colorbar
    shading flat
    if isfield(centroids,'admin0_name')
        title_string = sprintf('Evapotranspiration (mm/day), %s', centroids.admin0_name);
    else
        title_string = 'Evapotranspiration (mm/day)';
    end
    title(title_string)
    hold on
    set(gcf,'Color',[1 1 1]) % white figure background
end




