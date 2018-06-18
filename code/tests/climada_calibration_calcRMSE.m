function res = climada_calibration_calcRMSE(files,subS,snapS,rmv)

% Calculating of the RMSE for all shapefiles defined in files
% 
% MODULE:
%   flood
% NAME:
%   climada_calibration_calcRMSE
% PURPOSE:
%   
% CALLING SEQUENCE:
%   res = climada_calibration_calcRMSE(files)
% EXAMPLE:
%   climada_ls_multipleflow(lon,lat,elevation,25,1)
% INPUTS: 
%   files:   structure (files = dir(['path' '/*.shp']) which contains filepaths 
%            of shapefiles (fields: name and folder (generated by climada_calibration_ls, containing model-output)
%            with different flow path parameters sets (need to be defined in
%            filename)
%            shapefile need to include fields: length, area
%   subS:    original shapefile (generated by climada_calibration_subS2raster)
%            whith observed length and area fields. Need to have same
%            number of slides as model-output shape files
%       
% OPTIONAL INPUT PARAMETERS:
%   snapS:   shape file with snapped source areas. Need to be based on same
%            resolution as DEM of calibrated shape files. If snapS is provided Field
%            .max_srcslope is extracted and used to remove slide with a
%            maximum local source slope smaller than phi--> do not
%            propagate, stop immediately. also number of considered
%            observations is save in res.num_obs
%   rmv:     vector which includes info if slides need to be removed (1)
%            from the RMSE calculation
% OUTPUTS:
%   res:    
% MODIFICATION HISTORY:
% Thomas R�lli, thomasroelli@gmail.com, 20180523, init
% Thomas R�lli, thomasroelli@gmail.com, 20180604, option to remove slides
% Thomas R�lli, thomasroelli@gmail.com, 20180614, count number of
%  considered slides

%check arguments
if ~exist('files'), return; end
if ~exist('subS'), return; end
if ~exist('snapS'), snapS = []; end
if ~exist('rmv'), rmv = []; end

%init output
res_cali = repmat(struct('resol',[]),1,numel(files));

if ~isempty(snapS)
    max_slope = [snapS.max_srcslop];
end

for i=1:numel(files)
   fileparts = strsplit(files(i).name,'_');
   res_cali(i).resol = fileparts(1);
   %phi
   num_idx = regexp(fileparts{2},'\d');
   res_cali(i).phi = str2double(extractAfter(fileparts{2},num_idx(1)-1));
   %vmax
   num_idx = regexp(fileparts{3},'\d');
   res_cali(i).vmax = str2double(extractAfter(fileparts{3},num_idx(1)-1));
   %exp
   num_idx = regexp(fileparts{4},'\d');
   res_cali(i).exp = str2double(extractAfter(fileparts{4},num_idx(1)-1));
   %dH
   num_idx = regexp(fileparts{5},'\d');
   res_cali(i).dH = str2double(extractAfter(fileparts{5},num_idx(1)-1));
   %iT
   num_idx = regexp(fileparts{6},'\d');
   tmp_str = extractAfter(fileparts{6},num_idx(1)-1);
   res_cali(i).iT = str2double([tmp_str(1) '.' tmp_str(2:end)]);
   %perWTS
   num_idx = regexp([fileparts{7:end}],'\d');
   tmp_str = extractBetween([fileparts{7:end}],num_idx(1),num_idx(end));
   res_cali(i).perWT = tmp_str;
   
   %open corresponding caliS
   caliS = shaperead([files(i).folder filesep files(i).name]);
   
   %calculate RMSE for length
   obs = [subS.length];%observed lengths
   pred = [caliS.length];%predicted lengths
   
   %set removed to nan
   if ~isempty(snapS)
       obs(max_slope<res_cali(i).phi) = nan;
       pred(max_slope<res_cali(i).phi) = nan;
   end
   
   if ~isempty(rmv)
       obs(rmv) = nan;
       pred(rmv) = nan;
   end
   
   res_cali(i).rmse_lgt = sqrt(mean((obs-pred).^2,'omitnan'));
   
   %calculate RMSE for area
   obs = [subS.area];%observed area
   pred = [caliS.area];%predicted area
   res_cali(i).rmse_area = sqrt(mean((obs-pred).^2,'omitnan'));
   
   if ~isempty(snapS)
        res_cali(i).num_obs = numel(obs)-sum(max_slope<res_cali(i).phi);
   end
   
   %save source file
   res_cali(i).source = [files(i).folder filesep files(i).name];
end

res = res_cali;

end