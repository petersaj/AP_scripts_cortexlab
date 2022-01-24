function dlc = AP_load_dlc(fn)
% dlc = AP_load_dlc(fn)
% 
% Load and table-format DeepLabCut output

% Hard-coded options corresponding to output
opts = detectImportOptions(fn);
opts.VariableNamesLine = 2;
opts.VariableDescriptionsLine = 3;
opts.RowNamesColumn = 1;
dlc_table = readtable(fn,opts,'ReadVariableNames',true);

% Clear first column (it's frame index - can't exclude in readtable)
dlc_table(:,1) = [];

% Get unique variable names and turn 'descriptions' to structure
dlc_vars = dlc_table.Properties.VariableNames;
dlc_fields = unique(regexprep(dlc_table.Properties.VariableNames,'_\d+',''));

dlc = struct;
for curr_field = dlc_fields
    curr_field = curr_field{:};
    curr_subfields = find(contains(dlc_table.Properties.VariableNames,curr_field));
    for curr_subfield = curr_subfields
        dlc.(curr_field).(dlc_table.Properties.VariableDescriptions{curr_subfield}) = ...
            table2array(dlc_table(:,curr_subfield));
    end
end

