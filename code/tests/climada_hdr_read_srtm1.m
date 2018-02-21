function min_max_lon_lat = climada_hdr_read(tif_filename,SRTM1)
% climada
% MODULE:
%   etopo
% NAME:
%   climada_hdr_read
% PURPOSE:
%   Read hdr files that go together with tif files and get coordinates
%   information (upper left, lower left, upper right, lower right)
% CALLING SEQUENCE:
%    min_max_lon_lat = climada_hdr_read(tif_filename)
% EXAMPLE:
%   min_max_lon_lat = climada_hdr_read('...\etopo\data\srtm_19_10.tif')
% INPUTS:
%   tif_filename: a char or a cell with one or multiple tif-filenames (e.g.
%   ...etopo\data\srtm_19_10.tif)
% OUTPUTS:
%   min_max_lon_lat: a 1x4 matrix with [min_lon max_lon min_lat max_lat] data
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150723, init based on climada_90m_DEM by Gilles Stassen
%-

% init
min_max_lon_lat = [];

global climada_global
if ~climada_init_vars,return;end % init/import global variables


if ~exist('tif_filename', 'var'), tif_filename  = []; end
if ~exist('SRTM1', 'var'), SRTM1  = 0; end

% convert to cell if filename is a char
if ischar(tif_filename); tif_filename = {tif_filename}; end
    
n_tiles = numel(tif_filename);

% read hdr information, the lat/lon extremes of each tile (in any order)
extremes.lon = [];
extremes.lat = [];
for tile_i = 1:n_tiles
    filename = dir(tif_filename{tile_i});
    hdr_filename = strrep(tif_filename{tile_i},'tif','hdr');
    try 
        fid = fopen(hdr_filename);
        
        switch SRTM1
            case 0 %for SRTM3
                scale_check = 0;
                while ~feof(fid),
                    line = fgetl(fid);
                    if scale_check
                        scale = str2num(line);
                        dlon = scale(1); dlat = scale(2);
                        scale_check =0;
                    end
                    if strfind(line,'ModelPixelScaleTag')
                        scale_check = 1;
                    end

                    if strfind(line,'Upper Left')
                        loc_i = strfind(line, '(');
                        loc_f = strfind(line, ')');
                        UL = str2num(line(loc_i+1:loc_f-1));
                    end
                    if strfind(line,'Lower Left')
                        loc_i = strfind(line, '(');
                        loc_f = strfind(line, ')');
                        LL = str2num(line(loc_i+1:loc_f-1));
                    end
                    if strfind(line,'Upper Right')
                        loc_i = strfind(line, '(');
                        loc_f = strfind(line, ')');
                        UR = str2num(line(loc_i+1:loc_f-1));
                    end
                    if strfind(line,'Lower Right')
                        loc_i = strfind(line, '(');
                        loc_f = strfind(line, ')');
                        LR = str2num(line(loc_i+1:loc_f-1));
                    end
                end %while
                fclose(fid);
                extremes.lon = [extremes.lon UL(1) UR(1)];
                extremes.lat = [extremes.lat UL(2) LL(2)];
            case 1 %for SRTM1
                while ~feof(fid),
                    line = fgetl(fid);
                    if strfind(line,'NROWS')
                        split = strsplit(line);
                        nrows = str2num(char(split(2)));
                    end
                    if strfind(line,'NCOLS')
                        split = strsplit(line);
                        ncols = str2num(char(split(2)));
                    end
                    if strfind(line,'ULXMAP')
                        split = strsplit(line);
                        UL(1) = str2num(char(split(2)));
                        LL(1) = UL(1);
                    end
                    if strfind(line,'ULYMAP')
                        split = strsplit(line);
                        UL(2) = str2num(char(split(2)));
                        UR(2) = UL(2);
                    end
                    if strfind(line,'XDIM')
                        split = strsplit(line);
                        xdim = str2num(char(split(2)));
                        UR(1) = UL(1)+xdim*(ncols-1);
                        LR(1) = UR(1);
                    end
                     if strfind(line,'YDIM')
                        split = strsplit(line);
                        ydim = str2num(char(split(2)));
                        LL(2) = UL(2)-ydim*(nrows-1);
                        LR(2) = LL(2);
                     end
                end %while
                fclose(fid);
                extremes.lon = [extremes.lon UL(1) UR(1)];
                extremes.lat = [extremes.lat UL(2) LL(2)];
        end %end switch
        %break;
    end %try
end %n_tiles
 
% struct containing the lat/lon extremes of each tile (in any order)
min_max_lon_lat = [min(extremes.lon) max(extremes.lon) min(extremes.lat) max(extremes.lat)];

return

