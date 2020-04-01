% Generate revision figures for ctx-str paper 
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

% NOTE: these are just in order that I wrote them at the moment

%% Load in task data

% Load data

% (task)
% data_fn = 'trial_activity_choiceworld'; % Primary dataset
% data_fn = 'trial_activity_choiceworld_4strdepth'; % Depth-aligned striatum
% exclude_data = true;

% (passive)
data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
exclude_data = false;

% (unused at the moment)
% data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
% exclude_data = true;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

%%












