function fl_Climada_Aqueduct_compare(country_risk,report_filename)
% Compare flood damages with target model 
% MODULE:
%   flood
% NAME:
%   fl_Climada_Aqueduct_compare
% PURPOSE:
%   Compare flood damages calculated by Climada to the flood damage 
%   estimates of the Aqueduct Global Flood Analyzer developed by the World 
%   Resources Insitute (WRI), which can be accessed here:
%       http://www.wri.org/resources/maps/aqueduct-global-flood-analyzer
% PREREQUISITES:
%   Country data from the Aqueduct flood tool needs to be downloaded from 
%   here: 
%   http://www.wri.org/resources/data-sets/aqueduct-global-flood-risk-maps
%   (Click on the yellow button reading "3 DOWNLOADS", select "Data by 
%   Country" and place the unzipped files in
%       climada_module_flood/data/system)
% PREVIOUS STEP:
%   Generate a "country flood risk" struct that collects all the flood 
%   damage data for a list of countries. You can create such a struct 
%   using the script fl_countryrisk_generate
% CALLING SEQUENCE:
%   fl_Climada_Aqueduct_compare(country_risk,report_filename)
% EXAMPLE:
%   fl_Climada_Aqueduct_compare(country_risk_calc('Barbados')); % all in one
% INPUTS:
%   fl_country_risk: a structure with the results from 
%   fl_countryrisk_generate
% OPTIONAL INPUT PARAMETERS:
%   report_filename: the filename of the Excel file the comparison is 
%       written to. Prompted for if not given 
% OUTPUTS: produces an excel sheet
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150326, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('country_risk','var'),return;end
if ~exist('report_filename','var'),report_filename='';end

% PARAMETERS
%
% define the return periods we report damage for
% set =[] to report expected damage (ED) only
% (note that the ED will be denoted as return period=0 in the report)
DFC_return_periods=[5 10]; % [100 250]
%
% module data directory
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% Aqueduct Global Flood Analyzer shapefile
aqueduct_shapefile = [module_data_dir filesep 'system' filesep ...
    'aqueduct_global_flood_risk_data_by_country_20150304.shp'];
%-

% read or load aqueduct shapefile
[fP,fN,~] = fileparts(aqueduct_shapefile);
aqueduct_shapefile_mat = [fP filesep fN '.mat'];
if climada_check_matfile(aqueduct_shapefile,aqueduct_shapefile_mat)
    % .mat file exists, we just have to load it
    load(aqueduct_shapefile_mat);
elseif exist(aqueduct_shapefile,'file')
    % shapefile exists; read it
    fprintf('Reading Aqueduct Global Flood Analyzer data from %s\n',...
        aqueduct_shapefile)
    aqueduct_data = climada_shaperead(aqueduct_shapefile);
    save(aqueduct_shapefile_mat,'aqueduct_data');
else
    % Error message
    cprintf([206 50 50]/255,['Error: Missing Aqueduct shapefile %s.\n'...
        'Can''t proceed.\n'],aqueduct_shapefile);
    return
end

% prompt for report_filename if not given
if isempty(report_filename) % local GUI
    report_filename=[climada_global.data_dir filesep 'results' filesep ...
        'Climada_Aqueduct_cmp.xls'];
    [filename, pathname] = uiputfile(report_filename, 'Save report as:');
    if isequal(filename,0) || isequal(pathname,0)
        report_filename=''; % cancel
    else
        report_filename=fullfile(pathname,filename);
    end
end

n_entities=length(country_risk);

next_res=1;

% prepare header and print format
header_str='admin0(country);ISO3;admin1(state/province);admin1_code;value;peril;return_period;damage_Climada;damage_Aqueduct;Climada/Aqueduct';
format_str='%s;%s;%s;%s;%g;%s;%i;%g;%f\n';
header_str=strrep(header_str,';',climada_global.csv_delimiter);
header_str=strrep(header_str,';',climada_global.csv_delimiter);
if ~isempty(DFC_return_periods),DFC_exceedence_freq = 1./DFC_return_periods;end

% collect results
for entity_i=1:n_entities
    
    if isfield(country_risk(entity_i).res,'hazard') % country exposed
        
        n_hazards=length(country_risk(entity_i).res.hazard);
        for hazard_i=1:n_hazards
            
            res(next_res).country_name=country_risk(entity_i).res.country_name;
            res(next_res).country_ISO3=country_risk(entity_i).res.country_ISO3;
            
            if isfield(country_risk(entity_i).res.hazard(hazard_i),'admin1_name')
                res(next_res).admin1_name=country_risk(entity_i).res.hazard(hazard_i).admin1_name;
                res(next_res).admin1_code=country_risk(entity_i).res.hazard(hazard_i).admin1_code; % does also exist
            else
                res(next_res).admin1_name='';
                res(next_res).admin1_code='';
            end
            
            if ~isempty(country_risk(entity_i).res.hazard(hazard_i).EDS)
                
                res(next_res).Value   =country_risk(entity_i).res.hazard(hazard_i).EDS.Value;
                res(next_res).peril_ID=country_risk(entity_i).res.hazard(hazard_i).EDS.hazard.peril_ID;

                ED(next_res)=country_risk(entity_i).res.hazard(hazard_i).EDS.ED; % we need for sort later
                res(next_res).return_period=0;
                res(next_res).damage=ED(next_res);
                res(next_res).annotation_name=country_risk(entity_i).res.hazard(hazard_i).EDS.annotation_name;
                
                if ~isempty(DFC_return_periods)
                    % calculate a few points on DFC
                    [sorted_damage,exceedence_freq] = climada_damage_exceedence(...
                        country_risk(entity_i).res.hazard(hazard_i).EDS.damage,...
                        country_risk(entity_i).res.hazard(hazard_i).EDS.frequency);
                    
                    nonzero_pos     = find(exceedence_freq);
                    sorted_damage   = sorted_damage(nonzero_pos);
                    
                    RP_damage       = interp1(exceedence_freq,sorted_damage,DFC_exceedence_freq);
                    
                    for RP_i=1:length(RP_damage)
                        res(next_res+1)=res(next_res); % copy
                        ED(next_res+1)=ED(next_res); % copy
                        next_res=next_res+1;
                        res(next_res).return_period=DFC_return_periods(RP_i);
                        %%%% extract here the corresponding damage from
                        %%%% Aqueduct
                        country_position = structfind(aqueduct_data,...
                            'unit_name',upper(res(next_res).country_name));
                        field_string = sprintf('G10_bh_%s',...
                            num2str(DFC_return_periods(RP_i)));  
                        damage_aq(RP_i) = ...
                            getfield(aqueduct_data,{country_position},field_string);
                        res(next_res).damage_Aqueduct = damage_aq(RP_i);
                        res(next_res).damage=RP_damage(RP_i);
                        res(next_res).Climada_Aqueduct=RP_damage(RP_i)/damage_aq(RP_i);
                        
                        
                    end % RP_i
                end % ~isempty(DFC_return_periods)
                
            else
                ED(next_res)          =0;
                res(next_res).ED      =0;
                res(next_res).Value   =0; % possible to improve
                res(next_res).peril_ID=country_risk(entity_i).res.hazard(hazard_i).peril_ID;
                res(next_res).annotation_name='EMPTY';
                res(next_res).admin1_name='';
            end % ~isempty(EDS)
            next_res=next_res+1;
            
        end % hazard_i
        
    end % country exposed
    
end % entity_i

% print results table
ED_index=length(ED):-1:1; % unsorted
fprintf(header_str);

for ED_i=length(ED_index):-1:1 % to sort descending
    
    % print to stdout
    fprintf(format_str,...
        res(ED_index(ED_i)).country_name,...
        res(ED_index(ED_i)).country_ISO3,...
        res(ED_index(ED_i)).admin1_name,...
        res(ED_index(ED_i)).admin1_code,...
        res(ED_index(ED_i)).Value,...
        res(ED_index(ED_i)).peril_ID,...
        res(ED_index(ED_i)).return_period,...
        res(ED_index(ED_i)).damage,...
        res(ED_index(ED_i)).damage_Aqueduct,...
        res(ED_index(ED_i)).Climada_Aqueduct);
    
    % fill the table to write to the Excel file
    excel_data{length(ED_index)-ED_i+1,1}=res(ED_index(ED_i)).country_name;
    excel_data{length(ED_index)-ED_i+1,2}=res(ED_index(ED_i)).country_ISO3;
    excel_data{length(ED_index)-ED_i+1,3}=res(ED_index(ED_i)).admin1_name;
    excel_data{length(ED_index)-ED_i+1,4}=res(ED_index(ED_i)).admin1_code;
    excel_data{length(ED_index)-ED_i+1,5}=res(ED_index(ED_i)).Value;
    excel_data{length(ED_index)-ED_i+1,6}=res(ED_index(ED_i)).peril_ID;
    excel_data{length(ED_index)-ED_i+1,7}=res(ED_index(ED_i)).return_period;
    excel_data{length(ED_index)-ED_i+1,8}=res(ED_index(ED_i)).damage;
    excel_data{length(ED_index)-ED_i+1,9}=res(ED_index(ED_i)).damage_Aqueduct;
    excel_data{length(ED_index)-ED_i+1,10}=res(ED_index(ED_i)).Climada_Aqueduct;
end % ED_i

if ~isempty(report_filename)
    
    % try writing Excel file
    if climada_global.octave_mode
        STATUS=xlswrite(report_filename,...
            {'admin0(country)','ISO3','admin1(state/province)','admin1_code','value','peril','return_period','damage_Climada','damage_Aqueduct','Climada/Aqueduct'});
        MESSAGE='Octave';
    else
        [STATUS,MESSAGE]=xlswrite(report_filename,...
            {'admin0(country)','ISO3','admin1(state/province)','admin1_code','value','peril','return_period','damage_Climada','damage_Aqueduct','Climada/Aqueduct'});
    end
    
    if ~STATUS || strcmp(MESSAGE.identifier,'MATLAB:xlswrite:NoCOMServer') % xlswrite failed, write .csv instead
        %MESSAGE.message % for debugging
        %MESSAGE.identifier % for debugging
        [fP,fN]=fileparts(report_filename);
        report_filename=[fP filesep fN '.csv'];
        fid=fopen(report_filename,'w');
        fprintf(fid,header_str);
        for ED_i=length(ED_index):-1:1
            fprintf(fid,format_str,...
                res(ED_index(ED_i)).country_name,...
                res(ED_index(ED_i)).country_ISO3,...
                res(ED_index(ED_i)).admin1_name,...
                res(ED_index(ED_i)).admin1_code,...
                res(ED_index(ED_i)).Value,...
                res(ED_index(ED_i)).peril_ID,...
                res(ED_index(ED_i)).return_period,...
                res(ED_index(ED_i)).damage,...
                res(ED_index(ED_i)).damage_Aqueduct)

        end % ED_i
        fclose(fid);
        fprintf('Excel failed, .csv report written to %s\n',report_filename);
    else
        [STATUS,MESSAGE]=xlswrite(report_filename,excel_data, 1,'A2'); %A2 not to overwrite header
        fprintf('report written to %s\n',report_filename);
    end
end


return
