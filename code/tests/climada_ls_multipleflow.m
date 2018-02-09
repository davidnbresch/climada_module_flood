function mult_flow = climada_ls_multipleflow(centroids,hazard,exponent)

% Calculating of outflow proportion from each cell to its neighbours.
% Calculation according to Horton et al. (2013)
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   
% INPUTS: 
%   centroids: a climada centroids stucture (including topographical
%              information)
% OPTIONAL INPUT PARAMETERS:
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow 
%  
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion
% Thomas Rölli, thomasroelli@gmail.com, 20180208, implementation of flow
%   path

%remove afterwards; load centroids and hazard
%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat')
%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_distance.mat')

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('hazard', 'var'), hazard = []; end
if ~exist('exponent', 'var'), exponent = []; end

% PARAMETERS 
if isempty(exponent); exponent = 4; end

%get dimension of grid field from lon/lat coordinates
%and reshap needed vectors --> easier to handel in grid format than in
%vector; only possible for regular placed gridpoints
%also flip matrix up down to have high latitude in beginning of array
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = flipud(reshape(centroids.lon,n_lat,n_lon));
lat = flipud(reshape(centroids.lat,n_lat,n_lon));
elevation = flipud(reshape(centroids.elevation_m,n_lat,n_lon));

%calculate gradients from each cell to its 8 neighbours
gradients = climada_centroids_gradients(lon,lat,elevation);

% mutiplying by 1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 


%%% calculate sum of all outflow cells
gradients_sum = sum(gradients,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0


%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)(i= 1 to 8))^x

mult_flow = ((gradients).^exponent)./((gradients_sum).^exponent);
%mult_flow = (gradients.^exponent)/(gradients_sum.^exponent);


% Abflussanteil berechnen (da mit exponent die summe nicht mehr 1); ist 
% allerdings nicht teil multiflow algorithmus gemäss holton et al.; evt.
% wieder entfernen
sum_mult_flow = sum(mult_flow,3);
sum_mult_flow(sum_mult_flow==0) = 1; %prevent divison by 0
mult_flow =  mult_flow./sum_mult_flow;

%assessing flow path; starting from cells which are equal 1 in
%hazard.intensity
intensity = flipud(full(reshape(hazard.intensity(1,:),n_lat,n_lon)));
active_cells = intensity;
next_active_cells = intensity;
spread = double(intensity);

%shif matrix such that intensity is spread from center to neighbour-cells
%starting at 12 o'clock and proceeding clockwise
shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1]*-1;

for runs=1:100
active_cells = next_active_cells;
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
            for c=1:8 %iteration through shift matrix
                if mult_flow(j,i,c) > 0 
                    %disp([j i c]);
                    %prevent non-valible indices at boarder 
                    %if ((j+shift_matrix(c,1)>0) && (j+shift_matrix(c,1)<=n_lat))...
                            %&& ((i+shift_matrix(c,2)>0) && (i+shift_matrix(c,2)>0<=n_lon))
                        %spread value of center to neighbour cell only if
                        %new value greater than old
                        if spread(j,i)*mult_flow(j,i,c) > spread(j+shift_matrix(c,1),i+shift_matrix(c,2))
                            spread(j+shift_matrix(c,1),i+shift_matrix(c,2))=spread(j,i)*mult_flow(j,i,c);
                        end
                        next_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
                    %end
                end
            end
            active_cells(j,i)=0;
        end
    end
end
end

surf(lon,lat,elevation,spread,'LineStyle','none');
disp('hier');











end
