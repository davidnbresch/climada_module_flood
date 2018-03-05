function LCCS = climada_lc_get(rectangle)
% climada
% MODULE:
%   
% NAME:
%   climada_lc_get
% PURPOSE:
%   Read a landcover data set from the Climate Change Initiative (CCI) 
%   project by ESA (300m resolution, for 2015). NetCDF-file can be downloaded under:
%   http://maps.elie.ucl.ac.be/CCI/viewer/
%   The code extracts the needed data (defined by rectangle) from the 
%   downloaded netCDF file. netCDF shoud be placed in your climada_data folder
% CALLING SEQUENCE:
%   LCCS = climada_lc_get(rectangle)
% EXAMPLE:
%   LCCS = climada_lc_get([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05])
% INPUTS:
%   rectangle = a rectangle to define lon/lat box [minlon maxlon minlat
%   maxlat]
% OUTPUTS:
%  LCCS: 
%       .lc(i,j): landcover data set
%       .x(i,j): the longitude coordinates
%       .y(i,j): the latitude coordinates
%       .flag_values: land cover class code
%       .flag_meanings: description of land cover classes
%       .cmap: RGB colors for colormap, corresponding to lc-classes
%       .sourcefile: the source file (e.g. .../...300m-P1Y-2015-v2.0.7.nc')
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180301, init

LCCS=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%check arguments
if ~exist('rectangle', 'var'), return; end

%paths
filename = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.nc';
path = [climada_global.data_dir filesep 'landcover'];
fullpath = [path filesep filename];
url = 'http://maps.elie.ucl.ac.be/CCI/viewer/';

%check if netCDF file exists
if ~exist([fullpath])
    cprintf([1,0.5,0],'netCDF-file not found. Need to be donwloaded on %s and saved in %s.',url,path)
    return
end

%read lat and lon vectors
lat = ncread([fullpath],'lat');
lon = ncread([fullpath],'lon');

%get lat and lon indices which covers study area (defined by rectangle)
min_lat = find(lat<rectangle(3),1,'first');
max_lat = find(lat>rectangle(4),1,'last');
min_lon = find(lon<rectangle(1),1,'last');
max_lon = find(lon>rectangle(2),1,'first');

%get final lat lon vectors
lat = lat(max_lat:min_lat);
lon = lon(min_lon:max_lon);

%read land cover data within rectangle
lccs = ncread(fullpath,'lccs_class',[min_lon max_lat],[numel(lon) numel(lat)]);
%convert from signed int8 to unsigned
lccs(find(lccs<0)) = lccs(find(lccs<0))+256;
lccs = uint8(lccs);

%lat lon vectors (also lccs) need to be transformed to regualar climada form
lat = flipud(lat);
lon = lon';
lccs = flipud(lccs');

[lon lat] = meshgrid(lon,lat);

%flag values and meaning
all_flag_values = ncreadatt(fullpath,'lccs_class','flag_values');
all_flag_values = double(all_flag_values);
%convert from signed int8 to unsigned
all_flag_values(find(all_flag_values<0)) = all_flag_values(find(all_flag_values<0))+256;
all_flag_values = uint8(all_flag_values);
all_flag_meanings = strsplit(ncreadatt(fullpath,'lccs_class','flag_meanings'),' ');
all_cmap = [0 0 0;255 255 100;255 255 100;255 255 0;170 240 240;220 240 100;...
    200 200 100;0 100 0;0 160 0;0 160 0;170 200 0;0 60 0;0 60 0;0 80 0;...
    40 80 0;40 80 0;40 100 0;120 130 0;140 160 0;190 150 0;...
    150 100 0;120 75 0;150 100 0;255 180 50;255 204 204;255 235 175;...
    255 210 120;255 235 175;0 120 90;0 150 120;0 220 130;...
    195 20 0;255 245 215;220 220 220;255 245 215;0 70 200;255 255 255]./255;

%extract needed flag and color information
unique_values = unique(lccs)';
i = ismember(all_flag_values,unique_values);
flag_values = all_flag_values(i);
flag_meanings = all_flag_meanings(i);
cmap = all_cmap(i,:);


%write information into centroids
LCCS.flag_values = flag_values;
LCCS.flag_meanings = flag_meanings;
LCCS.cmap = cmap;
LCCS.lc = lccs;
LCCS.x = lon;
LCCS.y = lat;





end