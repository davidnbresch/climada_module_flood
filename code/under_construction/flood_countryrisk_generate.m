% script go generade "country flood risks"
% How to use: Specify the countries you want to calculate the flood risks 
% for in "country_list" (in the PARAMETERS section just below). 
% The script will then calculate a struct "fl_country_risk" that contains 
% the flood risk for all these countries and can be given as an input to
% the function fl_Climada_Aqueduct_compare, which will write an excel
% spreadsheet that compares the flood damages for given return
% periods calculated by CLIMADA to the corresponding estimates by the
% Aqueduct Global Flood Analyzer (see fl_Climada_Aqueduct_compare for
% details)
% NOTE: This is a very stupid script that just assumes that all the 
% entities and centroids needed in the calculation are present in their
% respective default folders; no error checking is done!
% Modification history: 
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150331, initial
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
%-

% initialize output
fl_country_risk = [];

% PARAMETERS
% country list
country_list={
    'Bangladesh'
    'Switzerland'
    };
check_plots = 0;
%
country_data_dir = climada_global.data_dir;

for country_i=1:length(country_list)
    [country_name,country_ISO3,shape_index] = climada_country_name(country_list{country_i});
    country_name_char = char(country_name); % as to create filenames etc., needs to be char
    
    fl_country_risk(country_i).res.country_name = country_name_char;
    fl_country_risk(country_i).res.country_ISO3 = country_ISO3;
    
    % define easy to read filenames
    centroids_file     = [country_data_dir filesep 'centroids' filesep country_ISO3 '_' strrep(country_name_char,' ','') '_centroids.mat'];
    entity_file        = [country_data_dir filesep 'entities' filesep country_ISO3 '_' strrep(country_name_char,' ','') '_entity.mat'];

    load(centroids_file)
    load(entity_file)
    
    %%%%% Generate flood hazard set here
    centroids = fl_centroids_prepare(centroids,15,'',0,0); %in case they have not been prepared yet
    save(centroids_file,'centroids');
    
    rf_hazard_save_file = [climada_global.data_dir filesep 'hazards' ...
        filesep country_ISO3 '_' country_name_char '_RF.mat'];
    if ~exist(rf_hazard_save_file, 'file')
        hazard_rf = climada_rf_hazard_set('', centroids, rf_hazard_save_file, check_plots);
    else
        load(rf_hazard_save_file)
    end
    fl_hazard_save_file = [climada_global.data_dir filesep 'hazards' ...
        filesep country_ISO3 '_' country_name_char '_FL.mat'];
    if ~exist(fl_hazard_save_file,'file')   
        hazard = climada_fl_hazard_set(hazard_rf,centroids,fl_hazard_save_file, check_plots);
    else
        load(fl_hazard_save_file)
    end
        
    fl_country_risk(country_i).res.hazard.peril_ID = 'FL';
    fl_country_risk(country_i).res.hazard.hazard_set_file = fl_hazard_save_file;
    fl_country_risk(country_i).res.hazard.entity_file = entity_file;
    
    % Damage calculation
    fl_country_risk(country_i).res.hazard.EDS = climada_EDS_calc(entity,hazard);
end