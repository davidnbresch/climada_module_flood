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
lccs = ncread('C:\Users\Simon Rölli\Desktop\climada\climada_data\landcover\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.nc',...
    'lccs_class',[min_lon max_lat],[numel(lon) numel(lat)]);

%get lat lon grid
%first lat lon vectors (also lccs) need to be transformed to regualar climada form
lat = flipud(lat);
lon = lon';
lccs = flipud(lccs');

[lon lat] = meshgrid(lon,lat);

%write information into centroids
flag_values = ncreadatt(fullpath,'lccs_class','flag_values');
flag_meanings = strsplit(ncreadatt(fullpath,'lccs_class','flag_meanings'),' ');

LCCS.flag_values = flag_values;
LCCS.flag_meanings = flag_meanings;
LCCS.lc = lccs;
LCCS.x = lon;
LCCS.y = lat;



disp('hier')





end