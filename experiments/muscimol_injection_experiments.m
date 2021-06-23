% Days and locations of muscimol injections/washouts
%
% (these are all hard-coded for each mouse here, there's no current way to
% associate something like this with a day's experiments)


%% Initialize structure

init_cell = cell(0);
muscimol = struct('animal',init_cell,'recordings',init_cell);

%% Corticostriatal mice

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP092';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP093';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP094';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP095';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP096';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP097';
muscimol(curr_animal_idx).recordings =  {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};

%% tetO mice

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP100';
muscimol(curr_animal_idx).recordings =  {...
    '2021-05-17','V1'; ...
    '2021-05-18','washout'; ...
    '2021-05-19','FRm'; ...
    '2021-05-20','washout'; ...
    '2021-05-21','DCS'; ...
    '2021-05-24','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP101';
muscimol(curr_animal_idx).recordings =  {...
    '2021-06-14','FRm'; ...
    '2021-06-15','washout'; ...
    '2021-06-16','FRm'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','V1'; ...
    '2021-06-21','DCS'; ...
    '2021-06-22','washout';};
    
curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP103';
muscimol(curr_animal_idx).recordings =  {...
    '2021-06-14','DCS'; ...
    '2021-06-15','washout'; ...
    '2021-06-16','V1'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','FRm'; ...
    '2021-06-21','washout'};

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP104';
muscimol(curr_animal_idx).recordings =  {...
    '2021-06-14','DCS'; ...
    '2021-06-15','washout'; ...
    '2021-06-16','FRm'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','V1'; ...
    '2021-06-21','washout'};



%% Sanity check
% (check something: all days exist? days account for all experiments before
% ephys and after retinotopy?)







