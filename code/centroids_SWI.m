function centroids=centroids_SWI(centroids, check_plots)
% Assign soil water indices (SWIs) to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_SWI_assign
% PURPOSE:
%   Determine soil water indices (SWIs) for given centroids by
%   reading a grayscale image of global annual mean SWI.
%   SWI describes the amount of soil moisture present out of how much could
%   potentially exist in the soil; i.e., the current water storage divided
%   by the available water storage. SWI ranges from 0 to 100% and is given
%   by:
%   SWI = (soil_moisture-wilting_point)/(field_capacity-wilting_point)
%   SWI therefore equals 0 at the wilting point and 1 at the field capacity
%   For more information on SWI see e.g. this paper by the Global Soil
%   Wetness Project:
%       https://www.jstage.jst.go.jp/article/jmsj1965/77/1B/77_1B_305/_pdf
%   Daily global SWI images derived can be downloaded here:
%   http://land.copernicus.eu/global/products/swi
%   The image needs to be placed in data/system of the climada module flood
%
%   NOTE: Since only daily data are provided under the link above, a
%       grayscale picture of the annual mean LAI needs to be produced
%       manually adding monthly pictures and taking the average.
% CALLING SEQUENCE:
%   centroids = centroids_SWI_assign(centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_SWI_assign(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, SWI is calculated even if the centroids
%   already have a field 'soil_water_index' (default is 0)
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'SWI' (soil water index), which contains the annual mean soil water index
%   (in %) for each centroid
% NOTE:
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150325, initial


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids),climada_centroids_load;      end
if ~exist('check_plots', 'var') || isempty(check_plots), check_plots = 0;       end

% locate the module's data folder
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% the file with the SWI data
SWI_file=[module_data_dir filesep 'system' filesep ...
    'SWI_global_mean_2013.png'];
%
%
% Missing data value, to be found in the documentation here:
% http://land.copernicus.eu/global/products/swi
missing_data_value = 255;
%
% north south extension of the SWI image (longitude is covered in the full
% range from -180 to 180)
min_South=-90; % degree
max_North= 90;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global SWI map)
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon), min(centroids.lat),...
    max(centroids.lon), max(centroids.lat)];

[fP,fN,~] = fileparts(SWI_file);
SWI_file_mat = [fP filesep fN '.mat'];
if climada_check_matfile(SWI_file,SWI_file_mat)
    % .mat file exists, we just have to load it
    load(SWI_file_mat);
elseif exist(SWI_file,'file')
    % GeoTIFF file exists; read it
    fprintf('reading soil water indices from %s\n',SWI_file)
    full_SWI_img = imread(SWI_file);
    save(SWI_file_mat,'full_SWI_img');
else
    % Error message
    cprintf([206 50 50]/255,['ERROR: missing image %s.\nCan''t proceed.\n'],SWI_file);
    return
end

full_SWI_img=full_SWI_img(end:-1:1,:); % switch for correct order in latitude (images are saved 'upside down')

% Crop the image to the centroids' bounding box (for speedup)
xx = 360*(1:size(full_SWI_img,2))/size(full_SWI_img,2)+(-180); % -180..180
yy = (max_North-min_South)*(1:size(full_SWI_img,1))/size(full_SWI_img,1)+min_South;
pos_x = find(xx>=bbox(1) & xx<=bbox(3));
pos_y = find(yy>=bbox(2) & yy<=bbox(4));
img = full_SWI_img(pos_y,pos_x);
xx = xx(pos_x);
yy = yy(pos_y);
[X,Y] = meshgrid(xx,yy); % construct regular grid

% Set the missing data to zero
img(img==missing_data_value) = 0;

% Determine the pixel closest to each centroid
n_centroids = length(centroids.centroid_ID);
X_1D = reshape(X,[1,numel(X)]);
Y_1D = reshape(Y,[1,numel(Y)]);
fprintf('assigning soil water indices to centroids...')
for centroid_i=1:n_centroids
    distances=climada_geo_distance(centroids.lon(centroid_i),...
        centroids.lat(centroid_i),X_1D,Y_1D);
    [~,min_dist_index] = min(distances);
    [img_row, img_col] = ind2sub(size(X),min_dist_index);
    centroids.SWI(centroid_i) = double(img(img_row, img_col));
end

% convert range of greyscale values into actual soil water indices
% ranging from 0 to 100
centroids.SWI = centroids.SWI./255*100;

fprintf(' done\n');

if check_plots
    % convert to double
    soil_wetis=double(img);
    % plot the image (kind of 'georeferenced')
    pcolor(X,Y,soil_wetis./255*100); colorbar
    shading flat
    if isfield(centroids,'admin0_name')
        title_string = sprintf('Soil water index (%), %s', ...
            centroids.admin0_name);
    else title_string = 'Soil water index (%)';
    end
    title(title_string)
    hold on
    set(gcf,'Color',[1 1 1]) % white figure background
end



