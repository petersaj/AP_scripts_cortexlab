%% Brain/striatum movie

presentation_folder = 'C:\Users\Andrew\Dropbox\CarandiniHarrisLab\210119_Oxford';

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

slice_spacing = 10;
structure_alpha = 1;

brain_fig = figure('Color','w','Position', [398,84,1157,898]);

brain_color = [1,0.6,0.6];
brain_av = av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1;
brain_3d = isosurface(permute(brain_av,[3,1,2]),0);
brain_3d_smoothed = smoothpatch(brain_3d,1,20);
brain_patch = patch('Vertices',brain_3d_smoothed.vertices*slice_spacing, ...
    'Faces',brain_3d_smoothed.faces, ...
    'FaceColor',brain_color,'EdgeColor','none','FaceAlpha',structure_alpha, ...
    'DiffuseStrength',0.6,'SpecularStrength',0,'AmbientStrength',0.5);

striatum_id = 574;
striatum_color = [0,0,0.8];
striatum_av = av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == striatum_id;
striatum_3d = isosurface(permute(striatum_av,[3,1,2]),0);
striatum_3d_smoothed = smoothpatch(striatum_3d,1,5);
striatum_patch = patch('Vertices',striatum_3d_smoothed.vertices*slice_spacing, ...
    'Faces',striatum_3d_smoothed.faces, ...
    'FaceColor',striatum_color,'EdgeColor','none','FaceAlpha',structure_alpha, ...
    'DiffuseStrength',0.6,'SpecularStrength',0,'AmbientStrength',0.5);

axis vis3d image off;
xlim([-100,1480])

set(gca,'Zdir','reverse')

% Set the properties for the mouse photo
view([-164,23]);
h = camlight('headlight');
set(h,'style','infinite');
lighting gouraud

camlight(h,'headlight');

% FOR FUTURE: rotate whole brain
view([-164,23]);
n_steps = 60;
[az_start,el_start] = view;
cam_steps = [linspace(az_start,az_start-360,n_steps); ...
    linspace(el_start,el_start,n_steps)];
movie_frames = struct('cdata',[],'colormap',[]);
for i = 1:length(cam_steps)
    view([cam_steps(1,i),cam_steps(2,i)]);
    camlight(h,'left');
    movie_frames(i) = getframe(brain_fig);
end
movie_file = [presentation_folder filesep 'brain_rotate.avi'];
writerObj = VideoWriter(movie_file);
writerObj.FrameRate = length(cam_steps)/3;
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);



%% Plot example single-experiment kernels

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

figure;imagesc(reshape(ephys_kernel_depth(5).k_px{1},461,[]))
colormap(brewermap([],'PRGn'));
caxis([-0.002,0.002]);
axis image off


%% Colored kernels over time

% Load data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time
framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_px_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_animal,'uni',false);


% Get mean kernels
ctx_str_k_px_cat = vertcat(ctx_str_k_px_animal{:}); 

ctx_str_k_px_task_mean = nanmean(cat(5,ctx_str_k_px_cat{:,1}),5);
ctx_str_k_px_notask_mean = nanmean(cat(5,ctx_str_k_px_cat{:,2}),5);

n_depths = size(ctx_str_k_px_task_mean,4);

% Plot colored kernels over time
figure;
for plot_k_px = 1:2
    switch plot_k_px
        case 1
            k_px = ctx_str_k_px_task_mean;
        case 2
            k_px = ctx_str_k_px_notask_mean;
    end
    
    % Plot center-of-mass color at select time points
    k_px_com = sum(k_px.*permute(1:n_depths,[1,3,4,2]),4)./sum(k_px,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    use_colormap = flipud(min(jet(255)-0.2,1));
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
            [1,n_depths])*size(use_colormap,1)),use_colormap);
    end
    
    plot_t = [-0.1:0.025:0.1];
    k_px_com_colored_t = ...
        permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
        [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
        size(k_px_com_colored,2),3),[2,3,4,1]);
    
    k_px_max = squeeze(max(k_px,[],4));
    k_px_max_t = ...
        permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
        plot_t),length(plot_t),size(k_px_max,1), ...
        size(k_px_max,2)),[2,3,1]);
    
    weight_max = double(prctile(k_px_max(:),99));
    subplot(2,1,plot_k_px);
    for t_idx = 1:length(plot_t)
        subplot(2,length(plot_t),(plot_k_px-1)*length(plot_t) + t_idx);
        p = image(k_px_com_colored_t(:,:,:,t_idx));
        set(p,'AlphaData', ...
            mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title([num2str(plot_t(t_idx)),' s']);
    end
    
end

%% Muscimol change index

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

fluor_kernelroi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{1},stimIDs{1},'uni',false));
fluor_kernelroi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{2},stimIDs{2},'uni',false));

%%% MUSCIMOL CHANGE OPTION 1
% Fit scaling factor (change) pre-post
n_exps = length(stimIDs{1});
str_muscimol_change = nan(n_exps,n_depths);
ctx_muscimol_change = nan(n_exps,n_depths);
for curr_depth = 1:n_depths
    for curr_exp = 1:n_exps
        str_muscimol_change(curr_exp,curr_depth) = ...
            mua_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            mua_postmuscimol_mean(curr_exp,:,curr_depth)';
        
        ctx_muscimol_change(curr_exp,curr_depth) = ...
            fluor_kernelroi_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            fluor_kernelroi_postmuscimol_mean(curr_exp,:,curr_depth)';
    end
end

%%% MUSCIMOL CHANGE OPTION 2
% Get difference in average stimulus response
t_stim = t >= 0 & t <= 0.2;
mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
str_muscimol_change = (mua_avg_postmuscimol-mua_avg_premuscimol);%./(abs(mua_avg_premuscimol)+abs(mua_avg_postmuscimol));

fluor_avg_premuscimol = permute(nanmean(fluor_kernelroi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postmuscimol = permute(nanmean(fluor_kernelroi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
ctx_muscimol_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol);%./(abs(fluor_avg_premuscimol)+abs(fluor_avg_postmuscimol));

% Plot time courses and change
figure;
p = gobjects(n_depths,3);
for plot_str = 1:n_depths
    
    p(plot_str,1) = subplot(n_depths,3,(plot_str-1)*n_depths+1); hold on;
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel(['Str ' num2str(plot_str)]);
    axis square
    
    p(plot_str,2) = subplot(n_depths,3,(plot_str-1)*n_depths+2); hold on;
    AP_errorfill(t,nanmean(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel('Cortex ROI');
    axis square
    
    p(plot_str,3) = subplot(n_depths,3,(plot_str-1)*n_depths+3);
    plot(ctx_muscimol_change(:,plot_str),str_muscimol_change(:,plot_str),'.k','MarkerSize',20)
    xlabel(['Cortex ROI (post-pre)']);
    ylabel(['Str ' num2str(plot_str) ' (post-pre)']);
    line(xlim,xlim,'color','k','linestyle','--');

end
linkaxes(p(:,1),'xy');
linkaxes(p(:,2),'xy');

% (str/ctx muscimol statistics)
disp('Striatum/cortex muscimol change correlation:')
for curr_depth = 1:n_depths
    [r,p] = corr(str_muscimol_change(:,curr_depth), ...
        ctx_muscimol_change(:,curr_depth), ...
        'rows','complete','type','pearson');
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(p) ' r = ' num2str(r)]);
end


disp('Striatum/cortex muscimol signrank:')
for curr_depth = 1:n_depths
    p = signrank(mua_avg_premuscimol(:,curr_depth),mua_avg_postmuscimol(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(p)]);
    p = signrank(fluor_avg_premuscimol(:,curr_depth),fluor_avg_postmuscimol(:,curr_depth));
    disp(['Ctx ' num2str(curr_depth) ' p = ' num2str(p)]);
end

%% Cell type training stim changes


data_fns = { ...
    'trial_activity_naive_sua', ...
    {'trial_activity_trainedPassive_sua.mat', ... % original
    'trial_activity_ctx_passive_sua.mat', ...     % + cortex ephys
    'trial_activity_muscimol_passive_sua.mat'}};  % muscimol group

stim_act_celltype_training = cell(size(data_fns));
stim_act_celltype_group = cell(size(data_fns));
stim_act_celltype_cat = cell(size(data_fns));

for curr_data = 1:length(data_fns)
    
    preload_vars = who;
    
    data_fn = data_fns{curr_data};
    
    % (turn on warnings)
    warning on;
    
    % Load data (saved as structure trial_data_all)
    trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    
    if ischar(data_fn)
        % Single dataset
        disp(['Loading ' data_fn '...']);
        load([trial_data_path filesep data_fn]);
        
    elseif iscell(data_fn)
        % Multiple datasets to merge (this is really dirty)
        
        preload_vars = who;
        
        % (load in all datasets)
        clear temp_data
        for curr_load_data = 1:length(data_fn)
            disp(['Loading ' data_fn{curr_load_data} '...']);
            temp_data{curr_load_data} = load([trial_data_path filesep data_fn{curr_load_data}],'trial_data_all');
            %         trial_data_all_split(curr_data) = temp_data.trial_data_all;
        end
        
        % (find shared fields - need to iterate intersect if > 2 datasets)
        split_fieldnames = cellfun(@(x) fieldnames(x.trial_data_all),temp_data,'uni',false);
        intersect_fieldnames = intersect(split_fieldnames{1},split_fieldnames{2});
        if length(data_fn) > 2
            for curr_intersect = 3:length(data_fn)
                intersect_fieldnames = intersect(intersect_fieldnames,split_fieldnames{curr_intersect});
            end
        end
        
        % (initialize combined structure)
        trial_data_all = cell2struct(cell(size(intersect_fieldnames)),intersect_fieldnames);
        data_animals = arrayfun(@(x) temp_data{x}.trial_data_all.animals,1:length(data_fn),'uni',false);
        trial_data_all.animals = horzcat(data_animals{:});
        
        % (concatenate experiment fields)
        experiment_fields = cellfun(@(curr_field) ...
            length(temp_data{1}.trial_data_all.(curr_field)) == ...
            length(temp_data{1}.trial_data_all.animals) && ...
            iscell(temp_data{1}.trial_data_all.(curr_field)) && ...
            any(cellfun(@(x) iscell(x),temp_data{1}.trial_data_all.(curr_field))),intersect_fieldnames);
        
        for curr_field = intersect_fieldnames(experiment_fields)'
            for curr_load_data = 1:length(data_fn)
                trial_data_all.(curr_field{:}) = ...
                    cat(1,trial_data_all.(curr_field{:}), ...
                    reshape(temp_data{curr_load_data}.trial_data_all.(curr_field{:}),[],1));
            end
        end
        
        % (grab non-experiment fields from the first dataset)
        % (NOTE: this assumes they're the same)
        for curr_field = intersect_fieldnames(~experiment_fields & ...
                ~strcmp(intersect_fieldnames,'animals'))'
            trial_data_all.(curr_field{:}) = ...
                temp_data{1}.trial_data_all.(curr_field{:});
        end
        
        % (if mixing protocols: ensure stimIDs convered into contrast*side)
        for curr_animal = 1:length(trial_data_all.trial_info_all)
            for curr_day = 1:length(trial_data_all.trial_info_all{curr_animal})
                
                curr_stim = trial_data_all.trial_info_all{curr_animal}{curr_day}.stimulus;
                if length(unique(curr_stim)) == 3 && all(unique(curr_stim) == [1;2;3])
                    % Stim [1,2,3] are coded stim IDs
                    curr_stim_recode = curr_stim;
                    % For both mpep and signals, ID is 1 = left, 2 = center, 3 = right
                    curr_stim_recode(curr_stim == 1) = -1;
                    curr_stim_recode(curr_stim == 2) = 0;
                    curr_stim_recode(curr_stim == 3) = 1;
                    
                    trial_data_all.trial_info_all{curr_animal}{curr_day}.stimulus = ...
                        curr_stim_recode;
                end
                
            end
        end
        
        % clear temp variables
        clearvars('-except',preload_vars{:},'data_fn','trial_data_all')
        
    else
        error('Unrecognized data_fn type')
    end
    
    
    % Julie's bugfixes
    % (combines 2 TAN categories, integrates old 'other' with UINs, puts short
    % waveforms that look like TANs into new 'other')
    removeMissingSpikesUnits = true;
    for iAnimal = 1:size(trial_data_all.goodUnits, 1)
        for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)
            
            % 1. if flag remove spikes missing, add these to the goodUnits
            if removeMissingSpikesUnits
                trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                    & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
            end
            
            % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs.
            if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
                for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
                    pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
                    pss_allcat2(iCell) = pss_allcat2temp(1);
                end
                allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
                largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40);
                
                fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
                uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
                tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
                msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
                shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
                
                trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
                trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
                trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
                trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
                trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
                
                clearvars pss_allcat2
            end
        end
    end
    
    
    % Find fields with experiment data (cell arrays with length animals)
    data_struct_fieldnames = fieldnames(trial_data_all);
    experiment_fields = cellfun(@(curr_field) ...
        length(trial_data_all.(curr_field)) == length(trial_data_all.animals) && ...
        iscell(trial_data_all.(curr_field)) && ...
        any(cellfun(@(x) iscell(x),trial_data_all.(curr_field))),data_struct_fieldnames);
    
    % Load pre-marked experiments to exclude and cut out bad ones
    if exist('exclude_data','var') && exclude_data
        exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
        exclude_fn = 'bhv_use_experiments';
        load([exclude_path filesep exclude_fn]);
        
        % Pull out used experiments for the animals loaded
        use_experiments_animals = ismember(trial_data_all.animals,{bhv_use_experiments.animals});
        use_experiments = {bhv_use_experiments(use_experiments_animals).use_experiments}';
        
        % (old, when multiple kinds of exclusions
        %     exclude_fn{1} = 'bhv_use_experiments';
        %     % exclude_fn{2} = 'expl_var_use_experiments';
        %     use_experiments_all = {};
        %     for curr_exclude = 1:length(exclude_fn)
        %         curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
        %         use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
        %     end
        %     use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
        %         1:size(use_experiments_all,2),'uni',false)';
        
        % Cut out bad experiments for any experiment data fields
        if ~isempty(use_experiments)
            for curr_field = data_struct_fieldnames(experiment_fields)'
                trial_data_all.(cell2mat(curr_field)) = cellfun(@(x,expts) ...
                    x(expts),trial_data_all.(cell2mat(curr_field)), ...
                    use_experiments,'uni',false);
            end
        end
    end
    
    % If any animals don't have any data - throw away
    nodata_animals = cellfun(@(x) isempty(x),trial_data_all.trial_info_all);
    trial_data_all.animals(nodata_animals) = [];
    for curr_field = data_struct_fieldnames(experiment_fields)'
        trial_data_all.(cell2mat(curr_field))(nodata_animals) = [];
    end
    
    % BUGFIX? If any days don't have data - throw away
    % (this isn't present in MUA dataset, don't know where this is from)
    for curr_animal = 1:length(trial_data_all.animals)
        nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
        if any(nodata_days)
            for curr_field = data_struct_fieldnames(experiment_fields)'
                trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
            end
        end
    end
    
    % Unpack data structure into workspace then throw away
    arrayfun(@(x) assignin('base',cell2mat(x),trial_data_all.(cell2mat(x))),data_struct_fieldnames);
    clear trial_data_all
    
    % Get sample rate and set "baseline" time
    sample_rate = 1/mean(diff(t));
    t_baseline = t < 0;
    
    % Get if this is a task dataset
    task_dataset = exist('outcome_all','var');
    
    % Concatenate trial info data
    trial_info_fields = fieldnames(trial_info_all{end}{end});
    trial_info_allcat = cell2struct(arrayfun(@(curr_field) ...
        cell2mat(cellfun(@(x) x(curr_field), ...
        cellfun(@struct2cell,vertcat(trial_info_all{:}),'uni',false))), ...
        1:length(trial_info_fields),'uni',false),trial_info_fields,2);
    
    % Concatenate wheel
    wheel_allcat = cell2mat(vertcat(wheel_all{:}));
    
    % % (movement from mousecam if exists: normalize and concatenate)
    % if exist('movement_all','var')
    %     movement_all_norm = cellfun(@(x) cellfun(@(x) x./nanstd(x(:)), ...
    %         x,'uni',false),movement_all,'uni',false);
    %     movement_allcat = cell2mat(vertcat(movement_all_norm{:}));
    % end
    
    %%% Get task/stim-relevant
    
    if task_dataset
        
        % Get trial information
        trial_stim_allcat = trial_info_allcat.stimulus; % old: diff(trial_info_allcat.stimulus,[],2);
        trial_choice_allcat = -(trial_info_allcat.response-1.5)*2;
        trial_outcome_allcat = trial_info_allcat.outcome;
        
        % Get reaction time and t index for movement onset
        move_t = trial_info_allcat.stim_to_move;
        [~,move_idx] = min(abs(move_t - t),[],2);
        
        % Get outcome time
        outcome_allcat = cell2mat(vertcat(outcome_all{:}));
        [~,outcome_idx] = max(any(outcome_allcat,3),[],2);
        outcome_t = t(outcome_idx)';
        
        % Get wheel velocity
        wheel_velocity_allcat = wheel_allcat;
        [max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
            (bsxfun(@times,wheel_velocity_allcat,trial_choice_allcat) > 0)),[],2);
        max_vel = max_speed.*trial_choice_allcat;
        wheel_velocity_allcat_move = wheel_velocity_allcat;
        
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(wheel_velocity_allcat,1)
            wheel_velocity_allcat_move(i,:,:) = circshift(wheel_velocity_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
        end
        
    elseif isfield(trial_info_allcat,'stimulus')
        
        if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
            % Use stim IDs
            trial_stim_allcat = trial_info_allcat.stimulus;
            % For both mpep and signals, ID is 1 = left, 2 = center, 3 = right
            trial_stim_allcat(trial_stim_allcat == 1) = -1;
            trial_stim_allcat(trial_stim_allcat == 2) = 0;
            trial_stim_allcat(trial_stim_allcat == 3) = 1;
        else
            % Passive choiceworld uses stimID = side*contrast
            trial_stim_allcat = trial_info_allcat.stimulus;
        end
        
    end
    
    
    % Choose split for data
    trials_allcat = size(wheel_allcat,1);
    trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Concatenate depths
    depth_allcat = cell2mat(cellfun(@(x) x(~isnan(x)), ...
        vertcat(allDepths{:}),'uni',false));
    
    % Concatenate groups and make logical vectors
    groups_allcat = cell2mat(cellfun(@(x) x(~isnan(x)), ...
        vertcat(allGroups{:}),'uni',false));
    
    % Concatenate good units
    good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
    good_units_allcat = cell2mat(vertcat(goodUnits{:})')';
    
    % (groups are cell type * depth: msn, fsi, tan, th, ??)
    n_aligned_depths = 3; % hardcoded: I think not stored
    n_celltypes = max(groups_allcat)./n_aligned_depths;
    if n_celltypes ~= 5
        error('Incorrect number of celltypes')
    end
    celltype_allcat = ceil(groups_allcat./n_aligned_depths);
    celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};
    
    domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
        n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);
    
    % Get maps for all cells
    % (load master U)
    load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
    % (use t = 0 kernels for all cells, hardcoded at the moment)
    use_k_frame = 5;
    ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
        vertcat(ctx_str_k_all{:})','uni',false));
    ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));
    
    
    % Deconvolve fluorescence and get kernel ROIs
    fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);
    
    n_vs = size(fluor_all{end}{end},3);
    kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
    load(kernel_roi_fn);
    fluor_kernelroi_deconv_exp = cellfun(@(x) ...
        permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
        permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
        vertcat(fluor_all{:}),'uni',false);
    
    
    % Get animals/days/cell number
    % (note this includes cortical cells, which are NaNs in allGroups)
    animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
    days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
    neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
        vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
    recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...
        num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));
    
    
    %%% Align by relative depth
    ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_align_fn = ['ephys_depth_align.mat'];
    load([ephys_align_path filesep ephys_align_fn]);
    
    % Pull out relative depths for each cell
    depth_aligned = cell(size(mua_all));
    for curr_animal = 1:length(mua_all)
        curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
            % (relative to end)
            depth_aligned{curr_animal}{curr_day} = ...
                allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})) - ...
                ephys_depth_align(curr_animal_idx).str_depth(curr_day,2);
            %         % (relative to start:end)
            %         depth_aligned{curr_animal}{curr_day} = ...
            %             (allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})) - ...
            %             ephys_depth_align(curr_animal_idx).str_depth(curr_day,1))./ ...
            %             diff(ephys_depth_align(curr_animal_idx).str_depth(curr_day,:));
        end
    end
    
    depth_aligned_allcat = cell2mat(horzcat(depth_aligned{:})');
    
    %%% Align by domain
    
    % Load domain alignment
    n_aligned_depths = 3;
    ephys_kernel_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_kernel_align_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
    load([ephys_kernel_align_path filesep ephys_kernel_align_fn]);
    
    % Loop through cells, pull out domains
    domain_aligned = cell(size(mua_all));
    for curr_animal = 1:length(mua_all)
        depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
            
            str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
            kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};
            
            % Get depth groups
            n_depths = round(diff(str_depth)/200);
            depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
            
            % Get domain depth boundaries
            kernel_match_boundary_idx = unique([1;find(diff(kernel_match) ~= 0)+1;n_depths+1]);
            kernel_match_depth_edges = depth_group_edges(kernel_match_boundary_idx);
            % (extend first and last forever - sometimes cells a little past
            % the border but classified as str probably because changed border
            % calculation slightly)
            kernel_match_depth_edges(1) = -Inf;
            kernel_match_depth_edges(end) = Inf;
            
            kernel_match_idx = kernel_match(kernel_match_boundary_idx(1:end-1));
            
            % Assign neurons to domains
            domain_aligned{curr_animal}{curr_day} = ...
                discretize(allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})), ...
                kernel_match_depth_edges,kernel_match_idx);
            
        end
    end
    
    domain_aligned_allcat = cell2mat(horzcat(domain_aligned{:})');
    
       
    %%%%% Average stim activity
    
    % Split data by experiment
    mua_exp = vertcat(mua_all{:});
    wheel_exp = vertcat(wheel_all{:});
    stim_exp = mat2cell(trial_stim_allcat,use_split,1);
    
    % Get quiescent trials
    wheel_thresh = 0.025;
    quiescent_trials = cellfun(@(wheel) ...
        ~any(abs(wheel(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2), ...
        wheel_exp,'uni',false);
    
    % Baseline subtract spikes
%     mua_exp_baselinesub = cellfun(@(act,quiescent) act - ...
%         nanmean(reshape(act(quiescent,t < 0,:),[],1,size(act,3)),1), ...
%         mua_exp,quiescent_trials,'uni',false);
%     %%%%%%%%%%%% TO NOT BASELINE-SUBTRACT
    mua_exp_baselinesub = mua_exp; warning('TESTING no baselinesub');
    
    % Get average activity with no wheel movement and 100%R stim
    stim_act = cellfun(@(stim,quiescent,act) ...
        permute(nanmean(act(stim == 1 & quiescent,:,:),1), ...
        [3,2,1]),stim_exp,quiescent_trials,mua_exp_baselinesub,'uni',false);
    
    % Get average activity for each celltype/domain/recording
    stim_act_cat = cell2mat(stim_act);
    [stim_act_celltype,group_char] = grpstats(stim_act_cat(good_units_allcat,:), ...
        [celltype_allcat(good_units_allcat),domain_aligned_allcat(good_units_allcat), ...
        recordings_allcat(good_units_allcat)],...
        {'mean','gname'});
    
    group = cellfun(@str2num,group_char);
    
    stim_act_celltype_training{curr_data} = stim_act_celltype;
    stim_act_celltype_group{curr_data} = group;
    
    % Store activity of all DMS cells
    stim_act_cat = vertcat(stim_act{:});
    
    stim_act_celltype_cat{curr_data}{1} = stim_act_cat( ...
        good_units_allcat & domain_aligned_allcat == 1 & celltype_allcat == 1,:);
    stim_act_celltype_cat{curr_data}{2} = stim_act_cat( ...
        good_units_allcat & domain_aligned_allcat == 1 & celltype_allcat == 2,:);
    stim_act_celltype_cat{curr_data}{3} = stim_act_cat( ...
        good_units_allcat & domain_aligned_allcat == 1 & celltype_allcat == 3,:);
    
    clearvars('-except',preload_vars{:},'t','n_aligned_depths','n_celltypes');
    
end

% Plot average stim timecourse before/after training
figure;
p = reshape(tight_subplot(n_aligned_depths,n_celltypes),[n_celltypes,n_aligned_depths])';
training_col = [0,0,0;1,0,0];
for curr_training = 1:2
    for curr_depth = 1:n_aligned_depths
        for curr_celltype = 1:n_celltypes
            
            curr_data = stim_act_celltype_training{curr_training}( ...
                all(stim_act_celltype_group{curr_training}(:,1:2) == ...
                [curr_celltype,curr_depth],2),:);
            
            axes(p(curr_depth,curr_celltype));
            AP_errorfill(t',nanmean(curr_data,1)', ...
                AP_sem(curr_data,1)',training_col(curr_training,:));
            
            xlim([-0.2,1])
            
        end
    end
end

for curr_celltype = 1:n_celltypes
    linkaxes(p(:,curr_celltype),'y');
end

% Get average stim activity before/after training
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);
stim_act_celltype_training_avg = cellfun(@(x) ...
    nansum(x(:,stim_avg_t_idx),2),stim_act_celltype_training,'uni',false);

stim_act_avg_grp = cellfun(@(grp,act) ...
    accumarray(grp,act,max(grp,[],1),@nanmean,NaN), ...
    stim_act_celltype_group,stim_act_celltype_training_avg,'uni',false);


% (statistics)
disp('Pre/post ranksum (MSN,FSI,TAN,UIN):');
for curr_depth = 1:3
    curr_p = nan(1,4);
    for curr_celltype = 1:4
        curr_p(curr_celltype) = ...
            ranksum(reshape(stim_act_avg_grp{1}(curr_celltype,curr_depth,:),[],1), ...
            reshape(stim_act_avg_grp{2}(curr_celltype,curr_depth,:),[],1));
    end
    disp(['Str ' num2str(curr_depth) ': ' num2str(curr_p)]);
end


% Plot change by celltype
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:4
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
        
        subplot(n_aligned_depths,4,(curr_depth-1)*4+curr_celltype);
        plotSpread(curr_data,'distributionColors',min(celltype_col(curr_celltype,:)+0.2,1))
        errorbar(cellfun(@nanmean,curr_data),cellfun(@AP_sem,curr_data), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
        set(gca,'XTickLabel',{'Untrained','Trained'}, ...
            'XTickLabelRotation',45)
        ylabel('Spikes/s');
    end
end



% Plot change by celltype
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

baseline_t = [-0.5,0];
stim_t = [0,0.15];

baseline_t_idx = t >= baseline_t(1) & t <= baseline_t(2);
stim_t_idx = t >= stim_t(1) & t <= stim_t(2);

% stim_act_celltype_cat_t_avg = cellfun(@(x) ....
%     (nanmean(x(:,stim_avg_t_idx),2) - nanmean(x(:,baseline_t_idx),2))./ ...
%     nanmean(x(:,baseline_t_idx),2), ...
%     vertcat(stim_act_celltype_cat{:}),'uni',false);
stim_act_celltype_cat_t_avg = cellfun(@(x) ....
    nanmean(x(:,stim_avg_t_idx),2), ...
    vertcat(stim_act_celltype_cat{:}),'uni',false);

figure;
for curr_celltype = 1:3
    
    curr_data = stim_act_celltype_cat_t_avg(:,curr_celltype);
    
    subplot(3,1,curr_celltype);
    plotSpread(curr_data,'distributionColors',min(celltype_col(curr_celltype,:)+0.2,1))
    set(gca,'XTickLabel',{'Untrained','Trained'}, ...
        'XTickLabelRotation',45)
    ylabel('Spikes/s');
end




%% Remove transparency, replace with equivalent opaque color

fig_ax = get(gcf,'Children');
for curr_ax = 1:length(fig_ax)
    ax_objs = get(fig_ax(curr_ax),'Children');
    for curr_obj = 1:length(ax_objs)
        try
            old_transparency = get(ax_objs(curr_obj),'FaceAlpha');
            set(ax_objs(curr_obj),'FaceColor', ...
                min(get(ax_objs(curr_obj),'FaceColor') + old_transparency,1));
            set(ax_objs(curr_obj),'FaceAlpha',1);           
        catch me
            continue
        end
    end
end
disp('Removed transparency');























