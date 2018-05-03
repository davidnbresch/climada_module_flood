function S = climada_ls_FPcalibration(orgS,orgLon,orgLat,orgStart,orgEnd,...
    lon,lat,elevation,field,save_name)

% Script to calibrate flow path parameters
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%     calculate: 1 if area lenght start end need to be calcualted
%     calculate2: 1 if 
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180502, init

global climada_global
if ~climada_init_vars, return; end

if ~exist('calculate') calculate=0; end
if ~exist('calculate2') calculate2=0; end

%init S
S = [];

n_lon = size(lon,2);
n_lat = size(lat,1);
start = zeros(n_lat,n_lon);
end_ = zeros(n_lat,n_lon);

%caculation of unitarea --> not considering slope
deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
unitarea = dx*dy;

%prepare grid to use in dsearchn()
grid = [lon(:) lat(:)];

%when get to coarser resolution --> 
%indices of start areas
idx_st = find(orgStart~=0); %indices of original start
points = [orgLon(idx_st) orgLat(idx_st)];
k_st = dsearchn(grid,points); 

idx_en = find(orgEnd~=0); %indices of original end
points = [orgLon(idx_en) orgLat(idx_en)]; %original points
k_en = dsearchn(grid,points); %get indices of nearest point to original points

maxfield='AREA_GIS';
for i=1:numel(orgS)
    S(i).(field) = orgS(i).(field); %transform field ID
    
    ik_st = find(orgStart(idx_st) == orgS(i).(field)); %find corresponding index of field ID
    ik_en = find(orgEnd(idx_en) == orgS(i).(field));
    %write start and end coordinates in Structure
    S(i).start = [lon(k_st(ik_st)) lat(k_st(ik_st))];
    S(i).end = [lon(k_en(ik_en)) lat(k_en(ik_en))];
    
    %calculate distance of slide in new raster (with slope)
    [Ist,Jst] = ind2sub(size(lon),k_st(ik_st));
    [Ien,Jen] = ind2sub(size(lon),k_en(ik_en));
    dI = abs(Ist-Ien)*dy;
    dJ = abs(Jst-Jen)*dx;
    dz = elevation(k_st(ik_st))-elevation(k_en(ik_en));
    try 
        lgt = double(sqrt(sqrt(dI^2+dJ^2)^2+dz^2));
    catch
        lgt = double(0);
    end
    %write transformed length in structure
    S(i).t_length = lgt;
    
end
start(k) = orgStart(idx);


%indices of end areas
idx = find(orgEnd~=0);
points = [orgLon(idx) orgLat(idx)];
k = dsearchn(grid,points);
end_(k) = orgEnd(idx);





%area of each raster cell --> considering slope
cell_area = climada_centroids_area(lon,lat,elevation,0);

%parameters for flowpath
dH = 0;
exponent = 25;
v_max = 4;
phi = 22;
delta_i = 0.0003;
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];

%needed data for flowpath
source = zeros(size(elevation));
elevation = fillsinks(elevation);
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH,1);
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%init of matrix to save flow paths ID and log(maximum intensity)
slides_ID = zeros(size(elevation));
slides_int = zeros(size(elevation));


%init structure to save results
S = [];

%dist = 
all_source = find(start~=0);

while sum(active_cells) > 0
    for i=1:n_lat
        for j=1:n_lon
            tem_sources = find(start~=0)
            %[I,J] = ind2sub(
        end
    end
    
end

for i = 1:numel(orgS)
    source(find(start==orgS(i).(field))) = 1; %find start of slide i and set as source
    spreaded = climada_ls_propagation(source,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt);
   
    %calculate area
    %without consideration of slope --> unitarea
    slide = find(spreaded~=0);
    raster_area = numel(slide)*unitarea;
    S(i).unitarea = raster_area;
    
    %with consideration of slope
    raster_area_slope = sum(cell_area(slide));
    S(i).slopearea = raster_area_slope;
    
    %add to 
    source(slide) = S(i).(field);
    forced_slides = forced_slides+source;
    source = source*0;
    i
end


end