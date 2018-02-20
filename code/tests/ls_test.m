function ls_test()
%function to test stuff... can be deleted afterwards

rect = [8-0.05 8+0.05 47-0.05 47+0.05];

srtm_min_lon_ndx = fix(rect(1));
srtm_max_lon_ndx = ceil(rect(2))-1;
srtm_min_lat_ndx = fix(rect(3));
srtm_max_lat_ndx = ceil(rect(4))-1;

[I,J] = meshgrid(srtm_min_lon_ndx:srtm_max_lon_ndx,srtm_min_lat_ndx:srtm_max_lat_ndx);

%assign N,S,W,E according to its coordinates --> included in nomenclatur of
%SRTM1 tiles filenames
NS = string(zeros(numel(I(:,1)),numel(I(1,:))));
WE = string(zeros(numel(I(:,1)),numel(I(1,:))));
NS(J>=0) = 'n';
NS(J<0) = 's';
WE(I>=0) = 'e';
WE(I<0) = 'w';

n_tiles = numel(I);

srtm_filename = cell(n_tiles,1);


for tile_i = 1:n_tiles
    srtm_filename{tile_i} = strcat(NS(tile_i),num2str(J(tile_i),'%02.0f'),'_',WE(tile_i),num2str(I(tile_i),'%03.0f'),'_1arc_v3.tif');
end

disp('hier')

end