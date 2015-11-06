function entity = climada_entity_crop(entity, bounding_box, entity_scale_factor)
% climada
% NAME:
%   climada_clip_centroids_entity
% PURPOSE:
%   Given an entity struct on country level (generated using 
%   climada_create_GDP_entity or climada_nightlight_entity) this function 
%   crops the entity struct to a bounding box and/or scales its resolution
%   (interpolating values).
% CALLING SEQUENCE:
%   entity = climada_entity_crop(entity, bounding_box, [entity_scale_factor])
% EXAMPLE:
% INPUTS:
%   entity:         The entity created by climada_create_GDP_entity or
%                   climada_nightlight_entity
%   bounding_box:   An array of size 4 bounding the region of interest,
%                   defined by [min_lon max_lon min_lat max_lat]
% OPTIONAL INPUT PARAMETERS:
%   cale_factor:    Can be a single number (scales by same factor along lat 
%                   and lon), or an array of size 2, such that the first
%                   index indicates the lon, and the second the lat scale 
%                   factors. Default value is set to 2 along both lat and lon.
% OUTPUTS:
%   entity:         Entity struct
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20141121
%   Gilles Stassen, gillesstassen@hotmail.com, 20141218 change variable
%                       names whole_world_borders.lon/lat -> shapes.X/Y
%   Gilles Stassen, gillesstassen@hotmail.com, 20141223 add check_country input arg.
%   Gilles Stassen, gillesstassen@hotmail.com, 20150220 centroids routine
%   removed, new function: climada_clip_centroids_entity -> climada_entity_crop
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('entity', 'var'),                 entity                  = [];     end
if ~exist('bounding_box', 'var'),           bounding_box            = [];     end
if ~exist('entity_scale_factor', 'var'),    entity_scale_factor     = 1;     end

% prompt for entity if not given
if isempty(entity),entity=climada_entity_load;end
if isempty(entity)
    fprintf('ERROR: Please provide an entity struct as input');
    return;
end

% If scale_factor is an array of two (x & y) values, scale x & y coords
% independently, otherwise, use same factor for both.

switch numel(entity_scale_factor)
    case 2
        e_s_f_x = entity_scale_factor(1);
        e_s_f_y = entity_scale_factor(2);
    case 1
        e_s_f_x = entity_scale_factor;
        e_s_f_y = entity_scale_factor;
end

tmp_lon = entity.assets.lon < bounding_box(2) & ...
    entity.assets.lon > bounding_box(1);
tmp_lat = entity.assets.lat < bounding_box(4) & ...
    entity.assets.lat > bounding_box(3);

% Trim entity structure to bounding box
flds = fieldnames(entity.assets);
no_ea_ori = numel(entity.assets.Value);
for i = 1 : numel(flds)
    if numel(entity.assets.(flds{i})) == no_ea_ori
        entity.assets.(flds{i}) = entity.assets.(flds{i})((tmp_lon & tmp_lat));
    end
end

entity = climada_entity_resolution_scale(entity,e_s_f_x,e_s_f_y);

