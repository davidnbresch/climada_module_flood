function [data, x, y] = climada_grid2array(data_grid, x_vector, y_vector)
% climada
% MODULE:
%   flood
% NAME:
%   climada_grid2array
% PURPOSE:
%   Restructure gridded data into singleton arrays. For consistency with
%   imagesc, the first index in data_grid corresponds to the x coords,
%   while the second index corresponds to the y coords. Used in
%   climada_read_srtm_DEM, read_APHRO_MA_V1101 and climada_rf_hazard_set
% CALLING SEQUENCE:
% EXAMPLE:
% INPUTS:
%   data_grid:  Gridded data  
% OPTIONAL INPUT PARAMETERS:
%   x_vector:   Monotonically in/decreasing vector defining the x values of
%               each of the columns of the data grid. If y_vector is left
%               empty, x_vector can also be a vector of size 4, defining
%               the corners of the grid as [min(x) max(x) min(y) max(y)]
%   y_vector:   Monotonically in/decreasing vector defining the x values of
%               each of the rows of the data grid.
% OUTPUTS:
%   data:       The original data structured as a singleton array 
%   x:          The x data structured as singleton array, same length as data array
%   y:          The y data structured as singleton array, same length as data array
% MODIFICATION HISTORY:
%   Gilles Stassen 20150107
%   Gilles Stassen 20150224 major cleanup
%-
data = [];      x = [];     y = []; %init

if ~exist('data_grid',      'var'),     return;             end
if ~exist('x_vector',       'var'),     x_vector = [];      end
if ~exist('y_vector',       'var'),     y_vector = [];      end

if numel(x_vector) == 4 && isempty(y_vector)
    reference_box = x_vector;
elseif ~isempty(y_vector) && ~isempty(x_vector)
    reference_box = [min(x_vector) max(x_vector) min(y_vector) max(y_vector)];
else
    reference_box = [];
end

switch ndims(data_grid)
    case 2
        data    = zeros(numel(data_grid),1);
        x   = zeros(numel(data_grid),1);
        y   = zeros(numel(data_grid),1);
        
        [size_x, size_y] = size(data_grid);

        if isempty(reference_box)
            reference_box = [1 size_x 1 size_y];
        end
        
        dx = (reference_box(2) - reference_box(1))/(size_x-1);
        dy = (reference_box(4) - reference_box(3))/(size_y-1);

        for j = 1 : size_y
            ndx = (j-1)*size_x;
            data(ndx + 1 : ndx + size_x)    = data_grid(:,j);
            y   (ndx + 1 : ndx + size_x,1)  = (size_y - j) * dx;
            x   (ndx + 1 : ndx + size_x,1)  = (0:size_x-1) .* dy;
        end

        x = x + reference_box(1);
        y = y + reference_box(3);
    case 3
        [size_x, size_y, size_t] = size(data_grid);
        
        if isempty(reference_box)
            reference_box = [1 size_x 1 size_y];
        end
        
        dx = (reference_box(2) - reference_box(1))/(size_x-1);
        dy = (reference_box(4) - reference_box(3))/(size_y-1);
        
        try
            data    = zeros(size_x*size_y,size_t);
        catch
            data    = sparse(size_x*size_y,size_t);
        end
        x   = zeros(size_x*size_y,1);
        y   = zeros(size_x*size_y,1);

        for j = 1 : size_y
            ndx = (j-1)*size_x;
            data(ndx + 1 : ndx + size_x,:)      = squeeze(data_grid(:,j,:));
            y   (ndx + 1 : ndx + size_x,1)      = (size_y - j) * dx;
            x   (ndx + 1 : ndx + size_x,1)      = (0:size_x-1) .* dy;
        end

        x = x + reference_box(1);
        y = y + reference_box(3);
    otherwise
        return;
end