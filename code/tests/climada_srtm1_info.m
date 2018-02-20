function [srtm_info,is_mat] = climada_srtm_info(centroidsORcountryORshapes,silent_mode,SRTM1)
% climada
% MODULE:
%   etopo
% NAME:
%   climada_srtm_info
% PURPOSE:
%   Get tile information for a specified country, centroids or a shape,
%   i.e. which tiles to download from srtm directory.
%
%   Data needs to be manually downloaded from from
%   http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp,
%   (or possibly also ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/,
%   http://droppr.org/srtm/v4.1/6_5x5_TIFs/)
%
%   See also https://dds.cr.usgs.gov/srtm/version2_1/SRTM3
%   https://dds.cr.usgs.gov/srtm/version2_1/Documentation/Quickstart.pdf
%   https://dds.cr.usgs.gov/srtm/version2_1/Documentation/SRTM_Topo.pdf
%   
%   Next call: climada_srtm_get
% CALLING SEQUENCE:
%   [srtm_info,is_mat] = climada_srtm_info(centroidsORcountryORshapes,silent_mode)
% EXAMPLE:
%   srtm_info = climada_srtm_info
%   srtm_info = climada_srtm_info('Netherlands')
%   srtm_info = climada_srtm_info([min_lon max_lon min_lat max_lat])
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   centroidsORcountryORshapes:  centroids with lat/lon, country_name or
%   shapes with lat/lon, it not specified prompted for a country_name
%   SRTM1:
% OUTPUTS:
%   srtm_tile: A struct with information which srtm_files to download from srtm
%        http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
%        .country_name:  country name of selected country, or region
%        .n_tiles:  number of tiles that are required for the selected region/country
%        .min_max_lon_lat: a 4x1 vector containing the srtm tile indices (e.g [18 19 10 10] for El Salvador)
%        .srtm_filename: a cell with filenames to be downloaded from srtm
%        .srtm_save_file: a char with the srtm filename (as .mat)
%   is_mat: =1 if matfile exists, =0 if it does not exist (=1 usually only
%       after climada_srtm_get has been called)
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150720, init based on climada_90m_DEM by Gilles Stassen
% david.bresch@gmail.com, 20160122, srtm folder moved, some fixes (removed hard-wired paths)
% david.bresch@gmail.com, 20160126, file info to stdout improved
% david.bresch@gmail.com, 20160513, error in latitude selection for southern hemisphere fixed, prompting for country name if invalid
% david.bresch@gmail.com, 20170323, issue with latitude for southern hemisphere patched, be careful, i.e. check output (as the tiles will be downloaded from www)
% Thomas Rölli, thomasroelli@gmail.com, 20180215, including SRTM1

srtm_info = [];
is_mat = 0;

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if exist(climada_global.map_border_file, 'file')
    load(climada_global.map_border_file) % contains shapes
end

if ~exist('centroidsORcountryORshapes', 'var'), centroidsORcountryORshapes = []; end
if ~exist('silent_mode', 'var'), silent_mode  = ''; end
if ~exist('SRTM1', 'var'), SRTM1  = 0; end
if isempty(silent_mode),silent_mode = 0;end

% PARAMETERS
%
% where the user can (manually) download the SRTM tiles
%srtm_address = 'ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/ ';
srtm_address = 'http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp';

% init
country_name = [];

% check for srtm folder in {climada_global.data_dir}/srtm (we can not store to the
% elevation_models module, as this creates troubles when updating code via github)
switch SRTM1
    case 0
        srtm_data_dir= [climada_global.data_dir filesep 'srtm'];
        if ~isdir(srtm_data_dir),mkdir(fileparts(srtm_data_dir),'srtm');end % create, should it not exist (no further checking)
    case 1
        srtm_data_dir= [climada_global.data_dir filesep 'srtm1'];
        if ~isdir(srtm_data_dir),mkdir(fileparts(srtm_data_dir),'srtm1');end % create, should it not exist (no further checking)
end

if ~isempty(centroidsORcountryORshapes)
    if isstruct(centroidsORcountryORshapes) && isfield(centroidsORcountryORshapes,'centroid_ID')
        % input is centroids
        centroids   = centroidsORcountryORshapes; clear centroidsORcountryORshapes
        if isfield(centroids,'countryname')
            country_name = centroids.country_name;
            [country_name, ~, occurrence] = unique(country_name);
            country_name = country_name(mode(occurrence));
            country_name = country_name{1};
        else
            country_name = [];
        end
        % define rectangle box
        rect = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)];
    elseif isstruct(centroidsORcountryORshapes) && (isfield(centroidsORcountryORshapes,'X') || isfield(centroidsORcountryORshapes,'lon'))
        % input is probably shapes
        shapes = centroidsORcountryORshapes; clear centroidsORcountryORshapes
        if isfield(shapes,'lon'), shapes.X =shapes.lon; shapes.Y = shapes.lat; end
        rect = [min(shapes.X) max(shapes.X) min(shapes.Y) max(shapes.Y)];
    elseif isnumeric(centroidsORcountryORshapes) && numel(centroidsORcountryORshapes) == 4
        % input is centroids_rect  
        centroids   = [];
        % define rectangle box
        rect        = centroidsORcountryORshapes; clear centroidsORcountryORshapes
    elseif ischar(centroidsORcountryORshapes)
        % input is country name
        centroids   = [];
        country_name= centroidsORcountryORshapes; clear centroidsORcountryORshapes
        [country_name,country_ISO3,shape_index] = climada_country_name(country_name);
        if isempty(shape_index)
            fprintf('ERROR: invalid country name \n')
            [country_name,country_ISO3,shape_index] = climada_country_name('Single');
        end
        rect = [min(shapes(shape_index).X) max(shapes(shape_index).X) min(shapes(shape_index).Y) max(shapes(shape_index).Y)];
    end
    
else
    % ask for country name
    centroids   = []; clear centroidsORcountry
    [country_name,country_ISO3,shape_index] = climada_country_name('Single');
    country_name = char(country_name);
    if isempty(country_name), return; end % error message already printed in climada_country_name
    rect = [min(shapes(shape_index).X) max(shapes(shape_index).X) min(shapes(shape_index).Y) max(shapes(shape_index).Y)];
end

% % country string for fprintf and fig title
% if exist('country_name','var') && ~isempty(country_name)
%     cntry_str = sprintf(' for %s',country_name'); 
% else
%     cntry_str = ''; 
% end

% conversion to srtm tile indices
switch SRTM1
    case 0 %SRTM3
        rect_buffer      = 0.5;     % set buffer to avoid missing data due to imperfect conversions below
        if rect(1) <-179.9,rect(1) = rect(1)+359.9; end
        srtm_min_lon_ndx = max(ceil(72 * (     rect(1)-rect_buffer +180)/(179.28+180.00)),1);
        srtm_max_lon_ndx = min(ceil(72 * (     rect(2)+rect_buffer +180)/(179.28+180.00)),72);
        if rect(3)<0 && rect(4)<0
            if ~silent_mode,fprintf('swapping min/max lat on Southern hemisphere\n');end
            rr=rect(4);
            rect(4)=rect(3);
            rect(3)=rr;
        end
        srtm_min_lat_ndx_ = max(ceil(24 * (60 - rect(3)-rect_buffer     )/( 60.00+ 57.83)),1);
        srtm_max_lat_ndx_ = min(ceil(24 * (60 - rect(4)+rect_buffer     )/( 60.00+ 57.83)),24);
        srtm_min_lat_ndx = min(srtm_min_lat_ndx_,srtm_max_lat_ndx_);
        srtm_max_lat_ndx = max(srtm_min_lat_ndx_,srtm_max_lat_ndx_);
        [I,J] = meshgrid(srtm_min_lon_ndx:srtm_max_lon_ndx,srtm_min_lat_ndx:srtm_max_lat_ndx);
        J=J(end:-1:1); % switch order as latitude is from -60..60, but tiles are numbered from North
        % number of tiles
        n_tiles = (1+abs(srtm_max_lon_ndx-srtm_min_lon_ndx))*(1+abs(srtm_max_lat_ndx-srtm_min_lat_ndx));
        % set srtm filenames
        srtm_filename = cell(n_tiles,1);
        for tile_i = 1:n_tiles
            srtm_filename{tile_i} = strcat('srtm_',num2str(I(tile_i),'%02.0f'),'_',num2str(J(tile_i),'%02.0f.zip'));
        end
    case 1 %SRTM1
        srtm_min_lon_ndx = fix(rect(1));
        srtm_max_lon_ndx = ceil(rect(2))-1;
        srtm_min_lat_ndx = fix(rect(3));
        srtm_max_lat_ndx = ceil(rect(4))-1;
        [I,J] = meshgrid(srtm_min_lon_ndx:srtm_max_lon_ndx,srtm_min_lat_ndx:srtm_max_lat_ndx);
        n_tiles = numel(I);
        
        %assign N,S,W,E according to its coordinates --> included in nomenclatur of
        %SRTM1 tiles filenames
        NS = char(zeros(numel(I(:,1)),numel(I(1,:))));
        WE = char(zeros(numel(I(:,1)),numel(I(1,:))));
        NS(J>=0) = 'n';
        NS(J<0) = 's';
        WE(I>=0) = 'e';
        WE(I<0) = 'w';
        for tile_i = 1:n_tiles
            srtm_filename{tile_i} = strcat(NS(tile_i),num2str(J(tile_i),'%02.0f'),...
                '_',WE(tile_i),num2str(I(tile_i),'%03.0f'),'_1arc_v3.tif');
        end
end
% create a specific filename that containts the min_max_lon_lats, to store the requested tile
% in order to speed-up subsequent calls
srtm_mat_filename = sprintf('srtm_%d_%d_%d_%d',srtm_min_lon_ndx,srtm_max_lon_ndx,srtm_min_lat_ndx,srtm_max_lat_ndx);
%%srtm_save_file = [climada_global.data_dir filesep 'results' filesep srtm_mat_filename '.mat'];
srtm_save_file = [srtm_data_dir filesep srtm_mat_filename '.mat'];
if exist(srtm_save_file,'file'), is_mat = 1; end

% set output structure
srtm_info.country_name = country_name;
srtm_info.n_tiles = n_tiles;
srtm_info.srtm_filename = srtm_filename;
srtm_info.min_max_lon_lat = [srtm_min_lon_ndx srtm_max_lon_ndx srtm_min_lat_ndx srtm_max_lat_ndx];
srtm_info.srtm_save_file = srtm_save_file;

% check if this file exists already

% create command window output
if ~silent_mode
    fprintf('SRTM DATA for %s (90m Digital Elevation Model)\n',country_name)
    if ~is_mat
        fprintf('Please download tif files from <a href="%s">%s</a> \n',srtm_address,srtm_address)
        fprintf(' Note: for Tyle X we use XX and for Tyle Y we use YY here\n')
        fprintf(' unzip the following files (srtm_XX_YY.zip)\n',srtm_data_dir,filesep,filesep)
        fprintf('  - %s\n', srtm_filename{:})
        fprintf(' and save as %s%ssrtm_XX_YY%ssrtm_XX_YY.tif\n',srtm_data_dir,filesep,filesep)
        % fprintf('\t save in %s \n\t and unzip\n', module_data_dir)
        % fprintf('<a href="%s">%s</a> \n', srtm_address, srtm_address)
    else
        fprintf('exists already as a mat-file (%s) \n',srtm_save_file)
    end
end

end % climada_srtm_info