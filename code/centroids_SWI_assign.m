function centroids=centroids_SWI_assign(centroids, check_plots, force_recalc)
% Assign soil wetness indices (SWIs) to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_SWI_assign
% PURPOSE:
%   Determine soil wetness indices (SWIs) for given centroids by
%   reading a grayscale image of global annual mean SWI. 
%   SWI describes the amount of soil moisture present out of how much could 
%   potentially exist in the soil; i.e., the current water storage divided 
%   by the available water storage. SWI ranges from 0 to 100%. 
%
%   Daily global SWI images derived can be downloaded here:
%   http://land.copernicus.eu/global/products/swi
%   The image needs to be placed in data/system of the climada module flood
%
%   NOTE: Since only daily data are provided under the link above, a
%       grayscale picture of the annual mean LAI needs to be produced 
%       manually adding monthly pictures and taking the average. 
% CALLING SEQUENCE:
%   centroids = centroids_SWI_assign(centroids, check_plots, force_recalc)
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
%   'soil_water_index', which contains the annual mean soil water index 
%   (in %) for each centroid
% NOTE: 
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150325, initial


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_SWI_assign, can''t proceed.\n']);
    return;
end
if ~exist('check_plots', 'var') || isempty(check_plots) 
    check_plots = 0;    end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;  end

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
% north south extension of the ET image (longitude is covered in the full 
% range from -180 to 180)
min_South=-90; % degree
max_North= 90;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global ET map)
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon), min(centroids.lat),...
    max(centroids.lon), max(centroids.lat)];
%-


% We only calculate ET values if the centroids do not come equipped
% with them or if force_recalc is set to 1
if ~isfield(centroids,'soil_wetness_index') || force_recalc
    [fP,fN,~] = fileparts(SWI_file);
    SWI_file_mat = [fP filesep fN '.mat'];
    if climada_check_matfile(SWI_file,SWI_file_mat)
        % .mat file exists, we just have to load it
        load(SWI_file_mat);
    elseif exist(SWI_file,'file')
        % GeoTIFF file exists; read it
        fprintf('Reading soil wetness indices from %s\n',SWI_file)
        full_SWI_img = imread(SWI_file);
        save(SWI_file_mat,'full_SWI_img');
    else
        % Error message
        cprintf([206 50 50]/255,['Error: Missing image %s.\n'...
            'Can''t proceed.\n',SWI_file]);
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
    fprintf('Assigning soil wetness indices to centroids...\n')
    for centroid_i=1:n_centroids
        distances=climada_geo_distance(centroids.lon(centroid_i),...
            centroids.lat(centroid_i),X_1D,Y_1D);
        [~,min_dist_index] = min(distances);
        [img_row, img_col] = ind2sub(size(X),min_dist_index);
        centroids.soil_wetness_index(centroid_i) = double(img(img_row, img_col));
    end
    
    % convert range of greyscale values into actual soil wetness indices
    % ranging from 0 to 100 
    centroids.soil_wetness_index = centroids.soil_wetness_index./255*100;
    
    cprintf([23 158 58]/255,['Successfully completed calculation of '...
        'soil wetness indices.\n']);
    
    if check_plots
        % convert to double 
        soil_wetis=double(img);
        % plot the image (kind of 'georeferenced')
        pcolor(X,Y,soil_wetis./255*100); colorbar
        shading flat
        if isfield(centroids,'admin0_name')
            title_string = sprintf('Soil wetness index (%), %s', ...
                centroids.admin0_name);
        else title_string = 'Soil wetness index (%)';
        end
        title(title_string)
        hold on
        set(gcf,'Color',[1 1 1]) % white figure background
    end

else
    % centroids already have a field 'evapotranspiration'
    cprintf([23 158 58]/255,['Skipped - centroids already have'...
        'soil wetness indices.\n']);
end % if isfield(centroids,'soil_wetness_index')
   
end % centroids_SWI_assign



