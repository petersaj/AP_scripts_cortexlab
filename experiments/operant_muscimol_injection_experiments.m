% Days and locations of muscimol injections/washouts
%
% (these are all hard-coded for each mouse here, there's no current way to
% associate something like this with a day's experiments)


%% Initialize structure

init_cell = cell(0);
muscimol = struct('animal',init_cell,'day',init_cell,'area',init_cell);

%% Corticostriatal mice

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP092';
curr_recordings = {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP093';
curr_recordings = {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP094';
curr_recordings = {...
    '2021-04-13','AM'; ...
    '2021-04-14','washout'; ...
    '2021-04-16','V1'; ...
    '2021-04-17','washout'; ...
    '2021-04-19','FRm'; ...
    '2021-04-20','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP095';
curr_recordings = {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP096';
curr_recordings = {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP097';
curr_recordings = {...
    '2021-04-27','V1'; ...
    '2021-04-28','washout'; ...
    '2021-04-29','FRm'; ...
    '2021-04-30','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);


%% tetO mice

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP100';
curr_recordings = {...
    '2021-05-17','V1'; ...
    '2021-05-18','washout'; ...
    '2021-05-19','FRm'; ...
    '2021-05-20','washout'; ...
    '2021-05-21','DCS'; ...
    '2021-05-24','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP101';
curr_recordings = {...
    '2021-06-14','FRm'; ... % wrong size pipette, redid next time
    '2021-06-15','washout'; ...
    '2021-06-16','FRm'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','V1'; ...
    '2021-06-21','DCS'; ...
    '2021-06-22','washout';};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);
    
curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP103';
curr_recordings = {...
    '2021-06-14','DCS'; ...
    '2021-06-15','washout'; ...
    '2021-06-16','V1'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','FRm'; ...
    '2021-06-21','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP104';
curr_recordings = {...
    '2021-06-14','DCS'; ...
    '2021-06-15','washout'; ...
    '2021-06-16','FRm'; ...
    '2021-06-17','washout'; ...
    '2021-06-18','V1'; ...
    '2021-06-21','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP105';
curr_recordings = {...
    '2021-07-07','FRm'; ...
    '2021-07-08','washout'; ...
    '2021-07-09','DCS'; ...
    '2021-07-12','washout'; ...
    '2021-07-13','V1'; ...
    '2021-07-14','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP106';
curr_recordings = {...
    '2021-07-07','V1'; ...
    '2021-07-08','washout'; ...
    '2021-07-09','DCS'; ...
    '2021-07-12','washout'; ...
    '2021-07-13','FRm'}; % really bad prep - unusable
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP107';
curr_recordings = {...
    '2022-01-26','washout'; ...
    '2022-01-31','washout'; ...
    '2022-02-01','V1'; ...
    '2022-02-02','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);

curr_animal_idx = length(muscimol)+1;
muscimol(curr_animal_idx).animal = 'AP108';
curr_recordings = {...
    '2022-01-26','washout'; ...
    '2022-01-31','washout'; ...
    '2022-02-01','V1'; ...
    '2022-02-02','washout'};
muscimol(curr_animal_idx).day =  curr_recordings(:,1);
muscimol(curr_animal_idx).area =  curr_recordings(:,2);


%% Sanity check
% (check something: all days exist? days account for all experiments before
% ephys and after retinotopy?)

disp('Checking that days + imaging exist for each mouse...')
for curr_animal = 1:length(muscimol)
    % Get all experiments for animal
    animal = muscimol(curr_animal).animal;    
    experiments = AP_find_experiments(animal);
    imaging_experiments = [experiments.imaging];
    
    % Make sure all the injection days exist
    check_days = ismember(muscimol(curr_animal).day,{experiments(imaging_experiments).day});
    if any(~check_days)
       error('%s %s: non-existent recording day',animal,muscimol(curr_animal).day{~check_days});
    end   
    AP_print_progress_fraction(curr_animal,length(muscimol));
end


%% Save

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = [data_path filesep 'muscimol.mat'];
save(save_fn,'muscimol');
disp(['Saved ' save_fn]);


%% Load muscimol experiments and get fluorescence stats

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

for curr_animal = 1:length(muscimol)
    
    animal = muscimol(curr_animal).animal;
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(muscimol(curr_animal).day)
        
        % Load muscimol imaging
        preload_vars = who;

        day = muscimol(curr_animal).day{curr_day};
        curr_experiments = AP_list_experiments(animal,day);
        
        use_protocol = 'AP_stimWheelRight'; 
        experiment = curr_experiments(find( ...
            strcmp({curr_experiments.protocol},use_protocol),1,'last')).experiment;

        % Load experiment
        AP_load_experiment;
        
        % Get downsampled pixel stats
        aUdf = AP_align_widefield(Udf,animal,day);
        
        U_downsample_factor = 5;
        aUd = imresize(aUdf,1/U_downsample_factor,'bilinear');
        
        px_d = AP_svdFrameReconstruct(aUd,fVdf);
        
        px_d_std = std(px_d,[],3);
        px_d_mad = mad(px_d,[],3);
        px_d_skew = skewness(px_d,[],3);
        
        % Upsample to regular resolution
        px_std = imresize(px_d_std,U_downsample_factor,'nearest');
        px_mad = imresize(px_d_mad,U_downsample_factor,'nearest');
        px_skew = imresize(px_d_skew,U_downsample_factor,'nearest');
        
        % Store in muscimol structure
        muscimol(curr_animal).px_avg{curr_day} = avg_im;
        muscimol(curr_animal).px_std{curr_day} = px_std;
        muscimol(curr_animal).px_mad{curr_day} = px_mad;
        muscimol(curr_animal).px_skew{curr_day} = px_skew;        
        
        % Prep next loop
        AP_print_progress_fraction(curr_day,length(muscimol(curr_animal).day));
        clearvars('-except',preload_vars{:});
        
    end
end

% Save the muscimol structure with pixel data
save(muscimol_fn,'muscimol');
disp(['Saved ' muscimol_fn]);


%% Plot pixel stats

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Plot all pixel stats
for curr_animal = 7:length(muscimol)
    curr_px = cat(3,muscimol(curr_animal).px_mad{:}); 
    curr_areas = muscimol(curr_animal).area;
    AP_image_scroll(curr_px,curr_area);axis image
    set(gcf,'name',muscimol(curr_animal).animal);
end

% Plot V1 muscimol change index
use_v1_muscimol = {'AP100','AP105','AP106','AP107','AP108'};
use_v1_muscimol_idx = find(contains({muscimol.animal},use_v1_muscimol));

v1_muscimol_change = cell(length(use_v1_muscimol_idx),1);
for curr_animal_idx = 1:length(use_v1_muscimol_idx)
    curr_animal = use_v1_muscimol_idx(curr_animal_idx);
    
    curr_areas = muscimol(curr_animal).area;
    
    curr_v1_muscimol_idx = find(strcmp(curr_areas,'V1'));
    curr_washout_idx = find(strcmp(curr_areas,'washout'));
    curr_v1_washout_idx = curr_washout_idx(find(curr_washout_idx > ...
        curr_v1_muscimol_idx,1));
    
    curr_v1_muscimol_px = muscimol(curr_animal).px_mad{curr_v1_muscimol_idx};
    curr_v1_washout_px = muscimol(curr_animal).px_mad{curr_v1_washout_idx};
    
    curr_v1_change = (curr_v1_muscimol_px - curr_v1_washout_px)./ ...
        (curr_v1_muscimol_px + curr_v1_washout_px);
    
    v1_muscimol_change{curr_animal_idx} = curr_v1_change; 
end

v1_muscimol_change_cat = cat(3,v1_muscimol_change{:});
AP_image_scroll(v1_muscimol_change_cat,{muscimol(use_v1_muscimol_idx).animal});
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
axis image;

figure;
v1_muscimol_change_avg = nanmean(v1_muscimol_change_cat,3);
h = imagesc(nanmean(v1_muscimol_change_avg,3));
set(h,'AlphaData',sum(isnan(v1_muscimol_change_cat),3) < 1);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
axis image off
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Average V1 muscimol change index');




