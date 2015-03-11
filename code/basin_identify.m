function basin_IDs = basin_identify(lon_pts,lat_pts,lon_polygons,lat_polygons,basin_names)
% Identify basins in which a set of points are located
% MODULE:
%   tbd
% NAME:
%   basin_identify
% PURPOSE: 
%   Identify basins in which a given set of points is located
%   For a set of longitude and latitude points (usually from centroids), 
%   determine if they are located in a given set of polygons, and if so,
%   assign a basin name (or ID) to them.
% CALLING SEQUENCE: 
%   basin_IDs = basin_identify(lon_pts,lat_pts,lon_polygons,lat_polygons,basin_names)
% INPUTS:
%   lon_pts [deg E]: longitude of input points [0,360)
%   lat_pts [deg N]: latitude of input points [-90,90]
%   lon_polygons [deg E]: cell array of longitude values of polygon vertices
%       (must be ordered either clockwise or counterclockwise)    
%   lat_polygons [deg N]: cell array of latitude values of polygon vertices
%       (must be ordered either clockwise or counterclockwise)  
%   basin_names: cell array of basin labels (one per cell in lon_polygons)
% OPTIONAL INPUTS:  
%
% OUTPUTS:
%   basin_IDs: cell array of basin labels at input points; 0 if
%       outside all input polygons
%
% NOTE 1: This function has been adapted from code written by
% Dan Chavas, CEE Dept, Princeton University
% See http://www.mathworks.com/matlabcentral/fileexchange/48360-basin-identify-m/content/basin_identify/basin_identify.m
%
% NOTE 2: The current implementation also allows for basin names that
% consist of strings (rather than integer basin IDs). However, since for
% the time being the function will be used mostly with basin numbers, 
% basin_IDs will be converted from cell to a numeric vector in the end,
% which is only possible if all contents of the cell array are of the same
% data type. Therefore, remove the last command ("cell2mat") in case you'd 
% like to attach the basin names to the input centroids as a cell
%
% MODIFICATION HISTORY
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150306, initial

% initialize output
basin_IDs = cell(size(lon_pts));

% loop over basins and check which centroids belong to them
for basin_i=1:length(basin_names)
   
    basin_str = basin_names{basin_i};
    lon_polys_temp = lon_polygons{basin_i};
    lat_polys_temp = lat_polygons{basin_i};
    lon_polys_temp = [lon_polys_temp lon_polys_temp(1)];   %to wrap around
    lat_polys_temp = [lat_polys_temp lat_polys_temp(1)];
        
    i_basin = inpolygon(lon_pts,lat_pts,lon_polys_temp,lat_polys_temp);
    if(sum(i_basin)>0)
        basin_IDs(i_basin) = {basin_str};
    end
    clear i_basin
    
end

basin_IDs(cellfun(@isempty,basin_IDs)) = {0};

% Convert the contents of the cell array basin_IDs into a single matrix.
basin_IDs = cell2mat(basin_IDs);


% if mapping toolbox is installed:
% if check_plot
%     h_pt = geoshow(lat_pts,lon_pts,'DisplayType','Point');
%     set(h_pt,'Marker','x','Color','g','MarkerSize',15)
%     textm(lat_pts,lon_pts,basin_pts)
%     title('Red X = test points')
%     
%     sprintf('Note: Lines on output plot do not align perfectly with meridians')
%     sprintf('so it may look like there are points that are in the wrong place. They are not.')
%     
%     %% Save plot %%%%%%%%%%%%%%%%%%
%     plot_filename = sprintf('basin_identify_example.pdf');
%     saveas(gcf,sprintf('%s',plot_filename),'pdf')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% end