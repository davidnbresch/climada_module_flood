function centroids=centroids_BD(centroids, check_plots)
% Assign bulk density (BD) values to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_BD
% PURPOSE:
%   Determine bulk density (BD) for given centroids by
%   reading a grayscale GeoTIFF of global annual mean BD.
%   The global average annual BD GeoTIFF image derived from http://soilgrids.org/
%
%   It needs to be placed in data/system of the climada module flood
%   The code automatically reads the GeoTIFF from that URL if it does not
%   exist in the data folder of the flood module (but this is VERY slow...).
%   bulk density values are in m
% CALLING SEQUENCE:
%   centroids = centroids_BD(centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_BD(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMBDERS:
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, BD is calculated even if the centroids
%   already have a field 'bulk density' (default is 0)
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'BD_kg_m3', which contains the annual mean bulk density
%   (in mm/day) for each centroid
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
% north south extension of the BD image (longitude is covered in the full
% range from -180 to 180)
min_lat=-90; % degree
max_lat= 90;
%
% prepare bounding box (for speedup; we only want to look at a section of
% the global BD map)
bbox = [min(centroids.lon), min(centroids.lat),max(centroids.lon), max(centroids.lat)];

format_str = '%s';

if check_plots
    fig = figure('color','w');
    colorbar
end

for sd_i = 1:6
    
    flds{sd_i} = sprintf('sd%i',sd_i);
    
    % the file with the BD data
    BD_file = [module_data_dir filesep 'system' filesep 'soil_data' filesep sprintf('BLD_sd%i_M_02_apr_2014.tif',sd_i)];
    % BD_file = sprintf(BD_file,sd_i);
    
    % URL where the BD file can be downloaded
    BD_file_URL  = ['ftp://soilgrids:soilgrids@ftp.soilgrids.org/data/recent/' sprintf('BLD_sd%i_M_02_apr_2014.tif.gz',sd_i)];
    % BD_file_URL = sprintf(BD_file_URL,sd_i);
    
    [fP,fN,~] = fileparts(BD_file);
    BD_file_mat = [fP filesep fN '.mat'];
    if exist(BD_file_mat,'file') %climada_check_matfile(BD_file,BD_file_mat)
        % .mat file exists, we just have to load it
        msg_str = sprintf('loading bulk density .mat files (%i/6)\n',sd_i);
        fprintf(format_str,msg_str)
        format_str = [repmat('\b',1,length(msg_str)) '%s'];
        load(BD_file_mat,flds{sd_i});
        if exist(flds{sd_i},'var')
            full_BD_img.(flds{sd_i}) = eval(flds{sd_i}); 
            clear -regexp ^sd\d{1}$ % Clear variables starting with "sd" and followed by 1 digit
        end
        if exist(BD_file,'file'), delete(BD_file); end % .tif file is large (>1GB!), so delete from disk
    elseif exist(BD_file,'file')
        % .tif file exists; read it
        msg_str = sprintf('reading bulk density data from %s\n',BD_file);
        fprintf(format_str,msg_str)
        format_str = [repmat('\b',1,length(msg_str)) '%s'];
        full_BD_img.(flds{sd_i}) = imread(BD_file);
        %fprintf('saving as .mat file...')
        save(BD_file_mat, '-struct', 'full_BD_img', flds{sd_i})
        %save(BD_file_mat,'full_BD_img');
        delete(BD_file)     % .tif file is large (>1GB!), so delete from disk
        %fprintf('done \n')
    else
        % Read file from the web
        msg_str = sprintf('reading bulk density data from web: %s\n',BD_file_URL);
        fprintf(format_str,msg_str)
        format_str = [repmat('\b',1,length(msg_str)) '%s'];
        gunzip(BD_file_URL,fP);
        full_BD_img.(flds{sd_i}) = imread(BD_file);
        %fprintf('saving as .mat file...')
        save(BD_file_mat, '-struct', 'full_BD_img', flds{sd_i})
        %save(BD_file_mat,'full_BD_img');
        delete(BD_file)     % .tif file is large (>1GB!), so delete from disk
        %fprintf('done \n')
    end
    
    %determine max density
    max_BD = max(max_BD,max(max(full_BD_img.(flds{sd_i}))));
    
    % full_BD_img.(flds{sd_i})=flipud(full_BD_img.(flds{sd_i})); % switch for correct order in latitude (images are saved 'upside down')

    % Crop the image to the centroids' bounding box (for speedup)
    lon             = 360*(1:size(full_BD_img.(flds{sd_i}),2))/...
        size(full_BD_img.(flds{sd_i}),2)+(-180); % -180..180
    lat             = (max_lat-min_lat)*(1:size(full_BD_img.(flds{sd_i}),1))/...
        size(full_BD_img.(flds{sd_i}),1)+min_lat;
    lon_crop_ndx    = lon>=bbox(1) & lon<=bbox(3);
    lat_crop_ndx    = lat>=bbox(2) & lat<=bbox(4);
    tmp_img         = uint16(full_BD_img.(flds{sd_i})); clear full_BD_img
    tmp_img         = flipud(tmp_img);
    img(:,:,sd_i)   = tmp_img(lat_crop_ndx,lon_crop_ndx);
    lon             = lon(lon_crop_ndx);
    lat             = lat(lat_crop_ndx);
    [LON,LAT]       = meshgrid(lon,lat); % construct regular grid
    
    % plot data at each depth interval
    if check_plots
        figure(fig)
        imagesc(lon,lat,img(:,:,1));
        set(gca,'YDir','Normal')
        depths = [2.5 10 22.5 45 80 150]; % hardwired from orig data source
        title(sprintf('Soil bulk density at %3.1f cm depth (kg/m^3)',depths(sd_i)));
    end
end
close fig
img                 = double(img);
mean_bulk_density   = mean(img,3); % take mean over depths

fprintf('assigning bulk density values to centroids...')
centroids.BD_kg_m3  = interp2(LON,LAT,mean_bulk_density,centroids.lon,centroids.lat,'cubic');
centroids.RD        = centroids.BD_kg_m3 ./max_BD; 
fprintf(' done\n');

if check_plots
    % convert to double (from uint16)
    bulk_density=griddata(centroids.lon,centroids.lat,centroids.BD_kg_m3,LON,LAT);%  double(img);
    % plot the image (kind of 'georeferenced')
    figure('color', 'w')
    hold on
    climada_plot_world_borders
    axis([bbox(1) bbox(3) bbox(2) bbox(4)])
    pcolor(LON,LAT,bulk_density);
    colorbar
    shading flat
    if isfield(centroids,'country_name')
        title_string = sprintf('bulk density (mm), %s', centroids.country_name{1});
    else
        title_string = 'bulk density (mm/day)';
    end
    title(title_string)
    hold off
end