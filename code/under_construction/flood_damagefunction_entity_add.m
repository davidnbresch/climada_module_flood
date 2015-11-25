% This script just replaces the damagefunctions of entities with the
% damagefunctions read from an entity template file. 
% In particular, it can be used to equip entities with the newly defined
% flood damagefunction in climada_module_flood\data\entity_template.xlsx

% path and filename of the entity template file to read the damagefunctions
% from
entity_template = climada_xlsread('no','C:\Users\S5ZQSN\Desktop\climada_meta\climada_modules\climada_module_flood\data\entity_template.xlsx','damagefunctions');
% regexp expression that describes the files (entities) where the
% replacement should be done
entity_file_regexp = 'C:\Users\S5ZQSN\Desktop\climada_meta\climada\data\entities\*_entity*';

% find the desired entity / entities
fP = fileparts(entity_file_regexp);
D_entity_mat = dir(entity_file_regexp);

% loop over entity files and replace damagefunctions
for file_i=1:length(D_entity_mat)
    
    entity_file_i = [fP filesep D_entity_mat(file_i).name];
    try
        load(entity_file_i)
        entity.damagefunctions = entity_template;
        fprintf('saving %s in %s (by %s)\n',D_entity_mat(file_i).name,fP,mfilename)
        save(entity_file_i,'entity')
        
    catch
        fprintf('skipped (invalid entity): %s\n',D_entity_mat(file_i).name);
        entity.assets=[]; % dummy
    end
    
end 