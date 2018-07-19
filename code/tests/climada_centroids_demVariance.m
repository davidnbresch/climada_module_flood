function vari = climada_centroids_demVariance(lon_coarse,lat_coarse,dem_coarse,lon_dense,lat_dense,dem_dense)

% 
% MODULE:
%  flood
% NAME:
%  climada_probSet_ls
% PURPOSE:
%  Calculates the varience of the slope degree in the coarse grid by
%  considering the variance of the slope degree of the grid cells in the dense grid 
%  which have a common nearest grid cell in the coarse grid.
%  Therefore, in a first step, for each cell in the dense grid, the nearest
%  grid cell within the coarse grid is searched (considering the lat/lon distance
%  and not y/x in meters). Next, the variance of slope degree is calculated
%  by including all slope degrees within the dense grid with the same
%  nearest coarse grid cell and the value assigned to the corresponding
%  coarse grid cell.
% CALLING SEQUENCE:
%   climada_centroids_demVariance(lon_coarse,lat_coarse,dem_coarse,lon2_dense,lat2_dense,dem2_dense)
% EXAMPLE:
%   
% INPUTS:
%   lon_,lat_,dem_coarse:    elevation and its coordinates in coarse grid
%                           (variance is assigned from dense grid)
%   
%   lon_,lat_,dem_dense:    elevation and its coordinates in dense grid
%                           (variance of dense grid is calculated)
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:   
%   vari:   Matrix with same dimension as coarse grid. Saves variance of
%           the slope degree when considering the variance within the dense
%           grid. Values at boarder are set to NaN.
% MODIFICATION HISTORY:
%  Thomas Rölli, thomasroelli@gmail.com, 20180716, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon_coarse') lon_coarse = []; end
if ~exist('lat_coarse') lat_coarse = []; end
if ~exist('dem_coarse') dem_coarse = []; end
if ~exist('lon_dense') lon_dense = []; end
if ~exist('lat_dense') lat_dense = []; end
if ~exist('dem_dense') dem_dense = []; end

if isempty(lon_coarse); return; end
if isempty(lat_coarse); return; end
if isempty(dem_coarse); return; end
if isempty(lon_dense); return; end
if isempty(lat_dense); return; end
if isempty(dem_dense); return; end

%get size of dense and coarse grid
size_coarse = size(lon_coarse);
size_dense = size(lon_dense);

%define coordinates points
coor_dense = [lon_dense(:) lat_dense(:)];
coor_coarse = [lon_coarse(:) lat_coarse(:)];

%caculate slope --> field to take the variance
slope = climada_centroids_slope(lon_dense,lat_dense,dem_dense);
field = slope;

%search for nearest grid point of each grid point in dense grid
idx_nearest = knnsearch(coor_coarse,coor_dense);

%calculate variance for dense grid points with common nearest coarse grid
%point
vari = accumarray(idx_nearest,field(:),[],@var);

vari = reshape(vari,size_coarse);

% set boarder variance to nan --> have unrealistic high values because
% knnsearch assignes to many cells for one coarse grid cell if two grids
% not perfectly fitting
% boarder cells
vari(size_coarse(1),:) = nan; %upper boarder
vari(1,:) = nan; %lower boarder
vari(:,1) = nan; %left boarder
vari(:,size_coarse(2)) = nan; %right boarder

%boxplot([vari_norm(:),vari2_norm(:)])
% figure
% surface(dem_coarse,vari,'LineStyle','none')

%slope_coarse = climada_centroids_slope(lon_coarse,lat_coarse,dem_coarse);

%sorte according to slope of coarse grid
%[~,idx_sor] = sort(slope_coarse(:));

%moving average
%m = movmean(vari(idx_sor),100,'omitnan');
%m = movmean(vari(idx_sor),50,'SamplePoints',slope_coarse(idx_sor));

%plot(slope_coarse(idx_sor),vari(idx_sor))

%plot(slope_coarse(:),vari(:),'.')
%histogram(slope_coarse(:),vari(:))




end
