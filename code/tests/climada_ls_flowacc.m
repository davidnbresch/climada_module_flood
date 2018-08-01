function FL_acc = climada_ls_flowacc(lon,lat,elevation,topoTB)
% Calculates flowaccumulation (upstream contributing area)
% MODULE:
%   flood
% NAME:
%	climada_ls_flowacc
% PURPOSE:
%   Calculate flowaccumulation for given raster, that is on a regular grid.
%   The function applies a multiple-flow-direction method which,
%   in contrast to the simple D8 method, allows the runoff to flow to
%   multiple neighbouring cells. The distribution of the flow is determined
%   based on the respective gradients between the central cell and its
%   neighbouring cells. For more information about the used algorithm see.
%       Freeman, T.G. (1991): Calculating catchment area with
%       divergent flow based on a regular grid;
%       doi:10.1016/0098-3004(91)90048-I
%       https://ac.els-cdn.com/009830049190048I/1-s2.0-009830049190048I-main.pdf?_tid=866d08a6-c54f-11e7-b499-00000aab0f01&acdnat=1510233280_58e175c46862affef6785a80491025fe
%   The code provides two versions of the calculation. One is including
%   functions of TopoToolbox. It is recommended to use this version because
%   it is very fast and handles flat areas. The other version was
%   implemented for CLIMADA but does not handle flat areas properly -->
%   handling of flat areas need to be implemented
% CALLING SEQUENCE:
%   flowacc = climada_ls_flowacc(lon,lat,z)
% EXAMPLE:
%   flowacc = climada_ls_TWI_calc(DEM,1) --> takes calculation with
%   TopoToolbox by default --> DEM is a GRIDobj
%   flowacc = climada_ls_TWI_calc(multflow,0) --> takes climada code, not
%   recommended --> multflow is a 8d multiple flow matrix
%   
% INPUTS:
%   lon/lat:    longitudinal/latitudinal coordinates in grid format
%               -->lat(i,j)
%   elevation:  the elevation according to coordinates in grid format -->
%               elevation(i,j). If the grid includes some NaNs it can lead
%               to problems
% OPTIONAL INPUT PARAMETERS:
%   topoTB:     1/0 if set to 1 the flowaccumulation is calculated with the
%               functions from TopoToolbox (recommended). 0 --> similar
%               approach but doesn't handle flat areas and is much slower.
% OUTPUTS:
%   centroids: centroids with additional fields:
%       .FL_acc:   assigns a number for flow accumulation to each grid cell
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180305, init
% Thomas Rölli, thomasroelli@gmail.com, 20180314, include functions form
%   TopoToolbox

%init flowaccumulation matrix
FL_acc = [];

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end
if ~exist('topoTB', 'var'), topoTB = []; end

% parameters
if isempty(lon), return; end
if isempty(lat), return; end
if isempty(elevation), return; end
if isempty(topoTB), topoTB = 1; end

%pre-processing: fill sinks, and NaNs --> in TopoToolbox
DEM = GRIDobj(lon,lat,elevation);
DEM.Z = deminpaint(DEM.Z);
DEM = fillsinks(DEM);

elevation = flipud(DEM.Z);

if topoTB
    DEM = GRIDobj(lon,lat,elevation);
    %uses geodesic for flat areas --> more sophisticated than routeflat
    FD = FLOWobj(DEM,'multi');
    %use routeflats for flat areas
    %M = flowdir(lon,lat,z,'type','multi','routeflats','geodesic');
    %FD = FLOWobj(M,'cellsize',DEM.cellsize,'size',DEM.size);
    A = flowacc(FD);
    FL_acc = flipud(A.Z);
    %with old TopoToolbox version 1.06
    %fl = flowdir(lon,lat,z,'routeflats','route');
    %a = flowacc(fl,[nlat nlon]);
    %FL_acc = a;
else 
    %multiple direction flow
    mult_flow = climada_ls_multipleflow(lon,lat,elevation,1.1,0,1);
    %version where multipleflow directions is adjusted in flat regions by
    %using functions in TopoToolbox (see routeflats)
    %verison where multiple flow directions is dereied by using dH to have outflow in flat
    %areas --> but leads to too high values in flat areas, and areas
    %donwstream to flat areas.
    %mult_flow = climada_ls_multipleflow(lon,lat,z,1.1,1,0);
    %to get from outflow to inflow --> shift matrix
    inflow_proportion = circshift(mult_flow,4,3);
    total_field = mult_flow(:,:,1)*0+1;
    field = mult_flow(:,:,1)*0+1;
    inflow_temp = mult_flow(:,:,1)*0;
    shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1];
    sumold = 0;
    
%     figure('units','normalized','outerposition',[0 0 1 1])
%     s = surface(total_field);
%     %view([30,30,50])
%     colorbar
    
    while sum(field(:),'omitnan')~=sumold && sum(field(:),'omitnan')~=0
        sumold = sum(field(:));
        for i=1:8
            inflow_temp(:,:,i) = circshift(field,shift_matrix(i,:)).*circshift(inflow_proportion(:,:,i),shift_matrix(i,:));
        end
        field = sum(inflow_temp,3);
        total_field = total_field + field;
        
%         s.CData = total_field;
%         pause(0.001)
%         sum(field(:))
    end
    
     FL_acc = total_field;
    
end
    




end