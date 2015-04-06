function centroids=centroids_ET_assign(centroids, check_plots, force_recalc)
% Assign Evapotranspiration (ET) values to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_ET_assign
% PURPOSE:
%   Determine Evapotranspiration (ET) for given centroids by
%   reading a grayscale GeoTIFF of global annual mean ET. 
%   The global average annual ET GeoTIFF image derived from 2010-2013 MODIS
%   data can be downloaded here:
%   ftp://ftp.ntsg.umt.edu/pub/MODIS/NTSG_Products/MOD16/MOD16A3.105_MERRAGMAO/Geotiff/MOD16A3_ET_2000_to_2013_mean.tif 
%   It needs to be placed in data/system of the climada module flood
%   The code automatically reads the GeoTIFF from that URL if it does not 
%   exist in the data folder of the flood module (but this is VERY slow...).
%   Evapotranspiration values are in mm/yr.
% CALLING SEQUENCE:
%   centroids = centroids_ET_assign(centroids, check_plots, force_recalc)
% EXAMPLE:
%   centroids = centroids_ET_assign(centroids)
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
%   'evapotranpiration', which contains the annual mean evapotranspiration 
%   (in mm/yr) for each centroid
% NOTE: 
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150320, initial


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_ET_assign, can''t proceed.\n']);
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
% the file with the ET data
ET_file=[module_data_dir filesep 'system' ...
    filesep 'MOD16A3_ET_2000_to_2013_mean.tif'];
%
% URL where the ET file can be downloaded
ET_file_URL = ['ftp://ftp.ntsg.umt.edu/pub/MODIS/NTSG_Products/'...
    'MOD16/MOD16A3.105_MERRAGMAO/Geotiff/MOD16A3_ET_2000_to_2013_mean.tif'];
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
min_South=-60; % degree
max_North= 80;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global ET map)
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon), min(centroids.lat),...
    max(centroids.lon), max(centroids.lat)];
%-


% We only calculate ET values if the centroids do not come equipped
% with them or if force_recalc is set to 1
if ~isfield(centroids,'evapotranspiration') || force_recalc
    [fP,fN,~] = fileparts(ET_file);
    ET_file_mat = [fP filesep fN '.mat'];
    if climada_check_matfile(ET_file,ET_file_mat)
        % .mat file exists, we just have to load it
        load(ET_file_mat);
    elseif exist(ET_file,'file')
        % GeoTIFF file exists; read it
        fprintf('Reading evapotranspiration data from %s\n',ET_file)
        full_ET_img = imread(ET_file);
        save(ET_file_mat,'full_ET_img');
    else
        % Read the GeoTIFF file from the web
        fprintf('Reading evapotranspiration data from web:\n\t%s\n',...
            ET_file_URL)
        full_ET_img = imread(ET_file_URL);
        save(ET_file_mat,'full_ET_img');
    end
    
    full_ET_img=full_ET_img(end:-1:1,:); % switch for correct order in latitude (images are saved 'upside down')
    
    % Set all 'special cases' defined in the parameters section to 0 (i.e.,
    % to black)
    full_ET_img(full_ET_img>65500) = 0;
    
    % Crop the image to the centroids' bounding box (for speedup)
    xx = 360*(1:size(full_ET_img,2))/size(full_ET_img,2)+(-180); % -180..180
    yy = (max_North-min_South)*(1:size(full_ET_img,1))/size(full_ET_img,1)+min_South;
    pos_x = find(xx>=bbox(1) & xx<=bbox(3));
    pos_y = find(yy>=bbox(2) & yy<=bbox(4));
    img = full_ET_img(pos_y,pos_x);
    xx = xx(pos_x);
    yy = yy(pos_y);
    [X,Y] = meshgrid(xx,yy); % construct regular grid
    
    % Determine the pixel closest to each centroid
    n_centroids = length(centroids.centroid_ID);
    X_1D = reshape(X,[1,numel(X)]);
    Y_1D = reshape(Y,[1,numel(Y)]);
    fprintf('Assigning evapotranspiration values to centroids...\n')
    for centroid_i=1:n_centroids
        distances=climada_geo_distance(centroids.lon(centroid_i),...
            centroids.lat(centroid_i),X_1D,Y_1D);
        [~,min_dist_index] = min(distances);
        [img_row, img_col] = ind2sub(size(X),min_dist_index);
        centroids.evapotranspiration(centroid_i) = double(img(img_row, img_col));
    end
    
    % convert range of greyscale values (uint16, i.e. values range from 0 
    % to 65535) into actual evapotranspiration (in mm/yr). 
    centroids.evapotranspiration = centroids.evapotranspiration.*0.1;
    
    cprintf([23 158 58]/255,['Successfully completed calculation of '...
        'evapotranspiration.\n']);
    
    if check_plots
        % convert to double (from uint16)
        evapo=double(img);
        % plot the image (kind of 'georeferenced')
        pcolor(X,Y,evapo.*0.1); colorbar
        shading flat
        if isfield(centroids,'admin0_name')
            title_string = sprintf('Evapotranspiration (mm/yr), %s', ...
                centroids.admin0_name);
        else title_string = 'Evapotranspiration (mm/yr)';
        end
        title(title_string)
        hold on
        set(gcf,'Color',[1 1 1]) % white figure background
    end

else
    % centroids already have a field 'evapotranspiration'
    cprintf([23 158 58]/255,['Skipped - centroids already have '...
        'evapotranspiration data.\n']);
end % if isfield(centroids,'evapotranspiration')
   
end % centroids_ET_assign



