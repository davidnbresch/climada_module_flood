function factor_of_safety = climada_ls_cell_failure(centroids, soil_moisture)
% compute the factor of safety
% MODULE:
%   flood
% NAME:
%   climada_ls_cell_failure
% PURPOSE:
%   Compute the factor of safety, given an appropriately prepared set of
%   centroids and a soil moisture matrix
% PREVIOUS STEP:
%   climada_ls_cell_failure
% CALLING SEQUENCE:
%   factor_of_safety = climada_ls_cell_failure(centroids,soil_moisture)
% EXAMPLE:
%   factor_of_safety = climada_ls_cell_failure(centroids,soil_moisture)
% INPUTS:
%   centroids:          centroids with the necessary fields (use
%                       climada_fl_centroids_prepare)
%   soil_moisutr:       matrix with soil moisture values for each event at
%                       each centroid (in mm)
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   factor_of_safety:   a 2d matrix the same size as the soil_moisture
%                       input.
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20150330 init
%   Gilles Stassen, 20150710, soil moisture as 2d matrix with one time dimension
%-

global climada_global
if ~climada_init_vars,  return; end;

if ~exist('centroids',      'var') ||  ~exist('soil_moisture',  'var')
    cprintf([1 0 0], 'ERROR: centroids and soil moisture are both required inputs')
    return;
end

n_events = size(soil_moisture,1);

BD_kg_m3    = repmat(centroids.BD_kg_m3,    [n_events 1]);
RD          = repmat(centroids.RD,          [n_events 1]);
WHC_mm      = repmat(centroids.WHC_mm,      [n_events 1]);
SD_m        = repmat(centroids.SD_m,        [n_events 1]);
slope_deg   = repmat(centroids.slope_deg,   [n_events 1]);

% gravity
g = 9.81;

% calculate soil moisture
water_sw        =   1000 * g;   % Calculate specific weight of water from density (1000 kg/m^3)
PD_kg_m3        =   2650;       % typical particle density in kg/m^3
porosity        =   1 - BD_kg_m3./PD_kg_m3; % wiki
%porosity        =   centroids.WHC_mm ./ centroids.SD_mm;    % Calculate porosity as ratio max soil storage to soil depth

% Calculate effective internal angle of friction (Meyerhoff 1956)
EIF             =   25.*RD + 25; % (emperical relationship)

% Calculate water table height (WTH) [temporal vector in metres]
WTH_m             =   SD_m - (WHC_mm - soil_moisture)./1000;  %soil_moisture .* (SD / FC);
% WTH             =   centroids.SD_m .* (centroids.WHC_mm ./ soil_moisture);  %soil_moisture .* (SD / FC);

% Calculate specific weight of solids from soil bulk density
solids_sw       =   BD_kg_m3 .* g;

% Calculate specific weight of moist soil
moist_soil_sw   =   solids_sw .* (1 - porosity);

% Calculate specific weight of saturated soil
sat_soil_sw     =   solids_sw .* (1 - porosity) + water_sw .* porosity;

% Cohesion
EC              =   1500 * repmat(mean(centroids.LAI,1),[n_events 1]); % soil cohesion from leaf area index (proxy)

% Calculate shear strength and stress
strength        =   EC + (moist_soil_sw.*(SD_m - WTH_m) + (sat_soil_sw - water_sw).*WTH_m).*tand(EIF).*(cosd(slope_deg)).^2;
stress          =   (moist_soil_sw.*(SD_m - WTH_m) + sat_soil_sw.*WTH_m) .* sind(slope_deg) .* cosd(slope_deg);

% Calculate failure probability from factor of safety
factor_of_safety             =   strength ./ stress;

if any(any(factor_of_safety <0))
    cprintf([1 0.5 0],'WARNING: points found with negative factor of safety, setting to NaN\n')
    factor_of_safety(factor_of_safety<0) = NaN;
end

