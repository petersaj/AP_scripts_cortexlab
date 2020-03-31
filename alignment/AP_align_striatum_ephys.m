function [str_depth,aligned_str_depth_group] = AP_align_striatum_ephys
% Align striatum ephys across experiments
%
% requires:
% str_align - 'depth' or 'kernel'
% 'depth' requires n_aligned_depths

%% Pull variables from base workspace
% (bad practice - but this used to be a script instead of a function and I
% wanted to not store all the intermediate variables because they're
% common)

animal = evalin('base','animal');
day = evalin('base','day');
spike_times = evalin('base','spike_times');
spike_templates = evalin('base','spike_templates');
spike_depths = evalin('base','spike_depths');
template_depths = evalin('base','template_depths');
verbose = evalin('base','verbose');
str_align = evalin('base','str_align');
n_aligned_depths = evalin('base','n_aligned_depths');

%% Get striatum boundaries

%%% Get correlation of MUA in sliding sindows
depth_corr_window = 100; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

max_depths = 3840; % (hardcode, sometimes kilosort2 drops channels)

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.01; % seconds
spike_binning_t_edges = nanmin(spike_times):spike_binning_t:nanmax(spike_times);

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

mua_corr = corrcoef(binned_spikes_depth');


%%% Estimate start and end depths of striatum

% % end of striatum: biggest (smoothed) drop in MUA correlation near end
% groups_back = 30;
% mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
% mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
% median_corr = medfilt1(nanmedian(mua_corr_end,2),3);
% [x,max_corr_drop] = min(diff(median_corr));
% str_end = depth_corr_bin_centers(end-groups_back+max_corr_drop);

% (new method)
% end of striatum: minimum correlation on dim 1 * dim 2
% (to look for the biggest dead space between correlated blocks)
groups_back = 20;
mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
mean_corr_dim1 = nanmean(mua_corr_end,2);
mean_corr_dim2 = nanmean(mua_corr_end,1);
mean_corr_mult = mean_corr_dim1.*mean_corr_dim2';
[~,mean_corr_mult_min_idx] = min(mean_corr_mult);
str_end = depth_corr_bin_centers(end-groups_back + mean_corr_mult_min_idx - 2); % err early: back up 2 (100 um)

% start of striatum: look for ventricle
% (by biggest gap between templates)
min_gap = 200;
sorted_template_depths = sort([0;template_depths]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
if max_gap > min_gap
    str_start = sorted_template_depths(max_gap_idx+1)-1;
else
    str_start = sorted_template_depths(2);
end

str_depth = [str_start,str_end];

%% Group striatal units by alignment

switch str_align
    
    case 'none'
        %%% Don't do anything
        
    case 'depth'
        %%% Align striatal recordings using saved alignment
        % (this requires an input for n_aligned_depths)
        % (NEW WAY: hardcode striatum length, doesn't need anything saved)
        if exist('n_aligned_depths','var')
            ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
            ephys_align_fn = ['ephys_depth_align.mat'];
            if exist([ephys_align_path filesep ephys_align_fn],'file')
                load([ephys_align_path filesep ephys_align_fn]);
                % Only use if this animal was included in the alignment
                curr_animal_idx = strcmp(animal,{ephys_depth_align.animal});
                if any(curr_animal_idx)
                    curr_day_idx = strcmp(day,ephys_depth_align(curr_animal_idx).day);
                    if any(curr_day_idx)
                        if verbose; disp('Aligning striatum by saved depths...'); end
                        
%                         % Get max striatum length and divide evenly
%                         all_str_lengths = diff(vertcat(ephys_depth_align(:).str_depth),[],2);
%                         max_str_length = max(all_str_lengths);
%                         global_str_depth_edges = linspace(0,max_str_length,n_aligned_depths+1);
                        
                        % ALTERNATE: set max striatum length, divide evenly
                        max_str_length = 3000; % Based on CCF trajectory
                        global_str_depth_edges = linspace(0,max_str_length,n_aligned_depths+1);
                        
                        % Set current striatum depth edges by global length
                        % assuming that the end is actually the end
                        str_depth_edges = sort(str_depth(2) - global_str_depth_edges);
                        str_depth_edges(1) = -Inf; % In case longer than max 
                        
                        % Get spike depths, setting all outside the striatum to NaN
                        str_spike_depths = spike_depths;
                        str_spike_depths(spike_depths < str_depth(1) | spike_depths > str_depth(2)) = NaN;
                        
                        % Group the striatal spike depths into the um-standardized bins
                        aligned_str_depth_group = discretize(str_spike_depths,str_depth_edges);
                        
                    end
                end
            end
        end
        
    case 'kernel'
        %%% Align by cortex-striatum kernel template
        ephys_kernel_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
        ephys_kernel_align_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];        
        if exist([ephys_kernel_align_path filesep ephys_kernel_align_fn],'file')
            load([ephys_kernel_align_path filesep ephys_kernel_align_fn]);
            % If alignment exists for this dataset, align
            curr_animal_idx = strcmp(animal,{ephys_kernel_align.animal});
            if any(curr_animal_idx)
                curr_day_idx = strcmp(day,ephys_kernel_align(curr_animal_idx).days);
                if any(curr_day_idx)
                    if verbose; disp('Aligning striatum by kernel alignment...'); end
                    % (use previously saved depth groups)
                    aligned_str_depth_group = ephys_kernel_align(curr_animal_idx).aligned_str_depth_group{curr_day_idx};
                    n_aligned_depths = ephys_kernel_align(curr_animal_idx).n_aligned_depths(curr_day_idx);
                end
            end
        end
        % Check that there's an aligned group for every spike
        if exist('aligned_str_depth_group','var') && length(spike_times) ~= length(aligned_str_depth_group)
            error('Not 1:1 raw and aligned spike depths')
        end
        
end

% If alignment isn't available, return empty
if ~exist('aligned_str_depth_group','var')
    aligned_str_depth_group = [];
end

%% Plot the aligned groups and MUA correlation

if verbose
    % Plot MUA correlation    
    figure;
    imagesc(depth_corr_bin_centers,depth_corr_bin_centers,mua_corr);
    axis tight equal;
    colormap(hot)
    line([str_depth(1),str_depth(1)],ylim,'color','b','linewidth',2);
    line([str_depth(2),str_depth(2)],ylim,'color','b','linewidth',2);
    line(xlim,[str_depth(1),str_depth(1)],'color','b','linewidth',2);
    line(xlim,[str_depth(2),str_depth(2)],'color','b','linewidth',2);
    xlabel('Probe depth (\mum)');
    ylabel('Probe depth (\mum)');
    title('MUA correlation: striatum location');
    drawnow;
end

if verbose && exist('aligned_str_depth_group','var') && ...
        ~isempty(aligned_str_depth_group)
    % Plot group depth
    [~,idx,~] = unique(spike_templates);
    template_aligned_depth = aligned_str_depth_group(idx)+1;
    template_aligned_depth(isnan(template_aligned_depth)) = 1;
    
    col = [0,0,1;copper(n_aligned_depths)];
    
    figure;plotSpread(template_depths,'distributionIdx', ...
        template_aligned_depth,'distributionColors',col(unique(template_aligned_depth),:));
    set(gca,'YDir','reverse');
    line(xlim,[str_depth(1),str_depth(1)]);
    line(xlim,[str_depth(2),str_depth(2)]);
    drawnow;
end

