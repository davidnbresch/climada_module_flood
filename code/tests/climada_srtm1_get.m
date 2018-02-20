function [SRTM,srtm_info] = climada_srtm_get(centroidsORcountryORshapes,check_plot,save_tile,verbose,SRTM1)
% climada
% MODULE:
%   elevation_models
% NAME:
%   climada_srtm_get
% PURPOSE:
%   Read the digital elevation model data from the files in an existing
%   srtm directory. Data can be downloaded from http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
%   or http://srtm.csi.cgiar.org/index.asp
%
%   The code tries to automatically get the tiles from
%   http://droppr.org/srtm/v4.1/6_5x5_TIFs/ (see SRTM_URL in PARAMETERS in
%   code). If this fails, the list of required tiles of the SRTM is printed
%   to stdout for the suer to manually retrieve
% CALLING SEQUENCE:
%   [SRTM,srtm_info] = climada_srtm_get(centroidsORcountryORshapes,check_plot,save_tile,verbose)
% EXAMPLE:
%   [SRTM,srtm_info] = climada_srtm_get('Switzerland',1,1,1)
%
%   SRTM = climada_srtm_get('El Salvador',1,0)
%   SRTM = climada_srtm_get([-89.15 -89.1 13.695 13.73],1,0) % las canas area, San Salvador
%   SRTM = climada_srtm_get([min_lon max_lon min_lat max_lat],1,0)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   centroidsORcountryORshapes:  can be a country_name 'El Salvador', or
%       coordinates of the rectangle box to get topography within [lonmin lonmax latmin latmax]
%       or centroids (as e.g. from climada_centroids_load)
%   check_plot: show a check plot, if =1, (default=0)
%   save_tile: =1: save the respective SRTM tiles, in order to speed up in
%       subsequent calls. =0 (default): do not save anything
%       Warning: the user is responsible for managing such tiles in the
%       ../data/results folder
%   verbose: =1 print info, =0 not (default)
%   SRTM1:
% OUTPUTS:
%   SRTM: a structure, with
%       .x(i,j): the longitude coordinates
%       .y(i,j): the latitude coordinates
%       .h(i,j): the elevation [m, negative below sea level]
%       sourcefile: the source file (e.g. .../srtm_18_10.tif)
%   If save_tile=1, the tile as returned is also saved to
%       ../data/results/etopo_*
% MODIFICATION HISTORY:
% muellele@gmail.com, 20150723, init based on climada_90m_DEM by Gilles Stassen and etopo_get by David Bresch
% david.bresch@gmail.com, 20160122, srtm folder moved, some fixes (removed hard-wired paths)
% david.bresch@gmail.com, 20160126, automatic retrieve implemented
% david.bresch@gmail.com, 20160126, single precision (half the memory need)
% david.bresch@gmail.com, 20160513, issue southern hemisphere solved
% david.bresch@gmail.com, 20160529, Cancel pressed works
% david.bresch@gmail.com, 20170323, issue with latitude for southern hemisphere patched, be careful, i.e. check output (as the tiles will be downloaded from www)
% Thomas Rölli, thomasroelli@gmail.com, 20180215, including SRTM1

SRTM=[];srtm_info=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('centroidsORcountryORshapes', 'var'), centroidsORcountryORshapes  = []; end
if ~exist('save_tile','var'), save_tile = [];  end
if ~exist('check_plot','var'), check_plot = 0; end
if ~exist('save_tile','var'), save_tile = 0; end
if ~exist('verbose','var'), verbose = 0; end
if ~exist('SRTM1', 'var'), SRTM1  = 0; end

% check for srtm folder in climada_data/srtm (we can not store to the
% elevation_models module, as this creates troubles when updating code via github)
switch SRTM1
    case 0
        srtm_data_dir= [climada_global.data_dir filesep 'srtm'];
        if ~isdir(srtm_data_dir),mkdir(fileparts(srtm_data_dir),'srtm');end % create, should it not exist (no further checking)
    case 1
        srtm_data_dir= [climada_global.data_dir filesep 'srtm1'];
        if ~isdir(srtm_data_dir),mkdir(fileparts(srtm_data_dir),'srtm1');end % create, should it not exist (no further checking)
end

% PARAMETERS
%
% the URL where to find the single tiles of the SRTM
switch SRTM1
    case 0
        SRTM_URL='http://droppr.org/srtm/v4.1/6_5x5_TIFs/';
    case 1
        SRTM_URL='https://earthexplorer.usgs.gov/';
end

if exist(climada_global.map_border_file,'file'),load(climada_global.map_border_file);end

% get srtm infos (filenames)
[srtm_info,is_mat] = climada_srtm1_info(centroidsORcountryORshapes,1,SRTM1);

if isempty(srtm_info),return;end % Cancel pressed
    
if is_mat
    % load from previously processed .mat file
    if verbose,fprintf('loading SRTM from %s\n',srtm_info.srtm_save_file);end
    load(srtm_info.srtm_save_file)
else
    % mat-file does not yet exist, so we create it
    % init SRTM structure, analog to ETOPO (see etopo_get.m)
    SRTM.x = [];
    SRTM.y = [];
    SRTM.h = [];
    SRTM.filename = srtm_info.srtm_save_file;
    SRTM.sourcefile = '';
    
    % check if the srtm files are downloaded and saved in climada data folder, subfolder srtm
    count = 0; %count number of SRTM1-files which need to be downloaded
    for tile_i = 1:srtm_info.n_tiles
        switch SRTM1
            case 0
                bare_filename=strrep(srtm_info.srtm_filename{tile_i},'.zip','');
                filename = [srtm_data_dir filesep bare_filename filesep bare_filename '.tif'];
            case 1
                bare_filename=strrep(srtm_info.srtm_filename{tile_i},'.tif','');
                filename = [srtm_data_dir filesep bare_filename '.tif'];
        end
        if ~exist(filename,'file')
            switch SRTM1
                case 0
                    % try to download and unzip
                    srtm_info.srtm_URL{tile_i} = [SRTM_URL srtm_info.srtm_filename{tile_i}];
                    try
                        fprintf('retrieving %s from %s (be patient) ... ',bare_filename,srtm_info.srtm_URL{tile_i});
                        unzip(srtm_info.srtm_URL{tile_i},fileparts(filename));
                        fprintf('done\n');

                        if ~exist(filename,'file')
                            fprintf('ERROR: retrieving %s worked, but unzip failed\n',srtm_info.srtm_filename{tile_i})
                            srtm_info = climada_srtm_info(centroidsORcountryORshapes,0);
                            return
                        end

                    catch
                        fprintf('ERROR: retrieving %s failed (likely sea tile, continuing)\n',srtm_info.srtm_filename{tile_i})
                        % new: continue (might be a sea tile with no elevation at all)
                        %srtm_info = climada_srtm_info(centroidsORcountryORshapes,0);
                        %return
                    end % try
                case 1
                    fprintf('%s not found!\n',srtm_info.srtm_filename{tile_i})
                    count = count+1;
            end
        end % ~exist(filename,'file')
        
    end % tile_i
    if SRTM1 && count>0
         cprintf([1,0.5,0],'%d file(s) not found (see above). Need to be donwloaded on %s and saved in %s.',count,SRTM_URL,srtm_data_dir)
        %return
    end
    
    % read srtm tif, including flip upside down
    % init
    if srtm_info.min_max_lon_lat(1)>srtm_info.min_max_lon_lat(2)
        srtm_info.min_max_lon_lat(1:2)=srtm_info.min_max_lon_lat(2:-1:1);
    end
    if srtm_info.min_max_lon_lat(3)>srtm_info.min_max_lon_lat(4)
        srtm_info.min_max_lon_lat(3:4)=srtm_info.min_max_lon_lat(4:-1:3);
    end
    
    [I,J] = meshgrid(srtm_info.min_max_lon_lat(1):srtm_info.min_max_lon_lat(2),...
        srtm_info.min_max_lon_lat(3):srtm_info.min_max_lon_lat(4));
    if SRTM1, k=I; I=J; J=k; end %switch because of different file terminology of SRTM1
    
    % srtm_raw_data.grid =
    for tile_i = 1:srtm_info.n_tiles
        switch SRTM1
            case 0
                bare_filename=strrep(srtm_info.srtm_filename{tile_i},'.zip','');
                filename = [srtm_data_dir filesep bare_filename filesep bare_filename '.tif'];
                [~, ~, fE] = fileparts(filename);
            case 1
                bare_filename=strrep(srtm_info.srtm_filename{tile_i},'.tif','');
                filename = [srtm_data_dir filesep filesep bare_filename '.tif'];
                [~, ~, fE] = fileparts(filename);
        end
        if exist(filename,'file')
            if strcmp(fE, '.tif') %|| strcmp(fE, '.tif.aux.xml')
                srtm_raw_data(I(tile_i),J(tile_i)).grid = flipud(imread(filename));
                SRTM.sourcefile{tile_i,1} = filename;
            end
        else
            srtm_raw_data(I(tile_i),J(tile_i)).grid = zeros(6000,6000); % patch empty
            SRTM.sourcefile{tile_i,1} = filename;
        end % exist(filename,'file')
    end % tile_i
    
    % Concatenate tiles to get one big tile which contains all the relavant
    % srtm info for one country, region or rectangle
    SRTM_grid = [];
    % loop over tiles in longitude direction (e.g. 18 and 19 for El Salvador)
    for i = srtm_info.min_max_lon_lat(1): srtm_info.min_max_lon_lat(2)
        SRTM_grid_j = [];
        % loopo over tiles in latitude direction (e.g. 10 for El Salvador)
        for j = srtm_info.min_max_lon_lat(3): srtm_info.min_max_lon_lat(4)
            SRTM_grid_j = [SRTM_grid_j ; srtm_raw_data(i,j).grid];
        end
        SRTM_grid = [SRTM_grid SRTM_grid_j];
    end
    clear srtm_raw_data SRTM_grid_j
    
    % overwrite no_data_value with a fill_data_value (either zero or the
    % real minimum data)
    % no_data value is the minimum of all values (i.e. -32768 for El Salvador)
    no_data_value = min(SRTM_grid(:));
    fill_data_value = nan; %0;
    %fill_data_value = min(SRTM_grid(SRTM_grid~= no_data_value));
    SRTM_grid(SRTM_grid == no_data_value) = fill_data_value;
    % convert to double
    SRTM_grid = double(SRTM_grid);
    
    % read lon/lat information from hdr
    min_max_lon_lat = climada_hdr_read(SRTM.sourcefile);
    
    % create meshgrid
    [x, y] = meshgrid(linspace(min_max_lon_lat(1),min_max_lon_lat(2),size(SRTM_grid,2)),...
        linspace(min_max_lon_lat(3),min_max_lon_lat(4),size(SRTM_grid,1)));
    
    % fill SRTM structure, analog to ETOPO (see etopo_get.m)
    % 20160126, single precision
    SRTM.x = single(x);
    SRTM.y = single(y);
    SRTM.h = single(SRTM_grid);
    
    if save_tile
        % save the SRTM for speedup in subsequent calls
        fprintf('saving SRTM as %s (delete to read new again)\n',SRTM.filename)
        save(SRTM.filename,'SRTM') % contains SRTM
    end
end %~is_mat

fprintf('SRTM area %f %f %f %f\n',min(min(SRTM.x)),max(max(SRTM.x)),min(min(SRTM.y)),max(max(SRTM.y)));

% cut out relevant area

if isstruct(centroidsORcountryORshapes)
    centroids_rect=[min(centroidsORcountryORshapes.lon) ...
        max(centroidsORcountryORshapes.lon) ...
        min(centroidsORcountryORshapes.lat) ...
        max(centroidsORcountryORshapes.lat)];
    clear centroidsORcountryORshapes
    centroidsORcountryORshapes=centroids_rect;
end % isstruct(centroidsORcountryORshapes)

if isnumeric(centroidsORcountryORshapes)
    if verbose,fprintf('restricting to [%2.3f %2.3f %2.3f %2.3f]... ',centroidsORcountryORshapes);end
    % we have a box of lon/lat
    % instead of inpoly, we can use a simple selector, as we have a
    % rectangle anyway
    %points = [SRTM.x(:) SRTM.y(:)];
    %polygon = centroidsORcountryORshapes([1 3; 2 4]);
    %is_selected = inpoly(points,polygon);
    is_selected = SRTM.x <= centroidsORcountryORshapes(2) & SRTM.x >= centroidsORcountryORshapes(1)...
        & SRTM.y <= centroidsORcountryORshapes(4) & SRTM.y >= centroidsORcountryORshapes(3);
    [I,J] = find(is_selected);
    min_I = min(I); max_I = max(I); min_J = min(J); max_J = max(J);
    SRTM.x = SRTM.x(min_I:max_I, min_J:max_J);
    SRTM.y = SRTM.y(min_I:max_I, min_J:max_J);
    SRTM.h = SRTM.h(min_I:max_I, min_J:max_J);
    if verbose,fprintf('done\n');end
end % isnumeric(centroidsORcountryORshapes)


if check_plot
    figure('Name','SRTM (90m digital elevation model)','Color',[1 1 1]);
    % create histogramm information to know relevant coloraxes
    %[no,xo] = hist(SRTM.h(:),0:500:7000);
    %no/sum(no)
    imagesc([min(SRTM.x(:)) max(SRTM.x(:))],[min(SRTM.y(:)) max(SRTM.y(:))],SRTM.h)
    set(gca,'YDir','normal')
    hold on
    colorbar;
    caxis;
    axis equal
    climada_plot_world_borders
    axis([min(SRTM.x(:)) max(SRTM.x(:)) min(SRTM.y(:)) max(SRTM.y(:))]);
    %climada_figure_axis_limits_equal_for_lat_lon([min(SRTM.x(:)) max(SRTM.x(:)) min(SRTM.y(:)) max(SRTM.y(:))])
    climada_figure_scale_add
end

% NOT properly implemented yet:
% smooth the DEM if desired
% if ~isempty(smooth) && any(smooth) && ~isnan(smooth)
%     if isscalar(smooth)
%         smooth_matrix = (1/smooth^2) .* ones(smooth);
%     elseif ismatrix(smooth)
%         smooth_matrix = smooth;
%     end
%     SRTM_grid = filter2(smooth_matrix,SRTM_grid);
% end

end % climada_srtm_get