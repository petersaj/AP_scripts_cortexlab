% AP_align_striatum_ephys
% Align striatum ephys across experiments
%
% requires:
% str_align - 'depth' or 'kernel'
% 'depth' requires n_aligned_depths

%% Get striatum boundaries

%%% Get correlation of MUA and LFP
n_corr_groups = 40;
max_depths = 3840; % (hardcode, sometimes kilosort2 drops channels)
depth_group_edges = linspace(0,max_depths,n_corr_groups+1);
depth_group = discretize(template_depths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges)-1;

spike_binning = 0.01; % seconds
corr_edges = nanmin(spike_times):spike_binning:nanmax(spike_times);
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

%%% Estimate start and end depths of striatum

% end of striatum: biggest (smoothed) drop in MUA correlation near end
groups_back = 15;
mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
median_corr = medfilt1(nanmedian(mua_corr_end,2),3);
[x,max_corr_drop] = min(diff(median_corr));
str_end = depth_group_centers(end-groups_back+max_corr_drop);

% start of striatum: look for ventricle (dropoff in templates)

% (by template density)
%     n_template_bins = 40;
%     size_template_bins = max(channel_positions(:,2))/n_template_bins;
%     template_density_bins = linspace(0,max(channel_positions(:,2)),n_template_bins);
%     template_density = histcounts(template_depths,template_density_bins);
%
%     str_end_bin = floor(str_end/size_template_bins);
%
%     n_bins_check = 3;
%     bins_conv = ones(1,n_bins_check)/n_bins_check;
%     template_gaps = conv(+(fliplr(template_density(1:str_end_bin)) < 2),bins_conv);
%
%     sorted_template_depths = sort([0;template_depths]);
%
%     if any(template_gaps)
%         str_gap_stop = length(template_gaps) - n_bins_check - find(template_gaps(n_bins_check:end),1);
%         str_start = sorted_template_depths(find(sorted_template_depths > template_density_bins(str_gap_stop),1)) - 1;
%     else
%         str_start = sorted_template_depths(2);
%     end

% (by biggest gap)
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
                        
                        % Get the maximum striatal length, split evenly
                        all_str_lengths = diff(vertcat(ephys_depth_align(:).str_depth),[],2);
                        max_str_length = max(all_str_lengths);
                        global_str_depth_edges = linspace(0,max_str_length,n_aligned_depths+1);
                        
                        % Set current striatum depth edges by global length
                        % assuming that the end is actually the end
                        str_depth_edges = sort(str_depth(2) - global_str_depth_edges);
                        
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
        % If the number of spike depths doesn't match depth groups, error
        if exist('aligned_str_depth_group','var') && length(spike_depths) ~= length(aligned_str_depth_group)
            error('Not 1:1 raw and aligned spike depths')
        end
        
end

%% Plot the aligned groups

if verbose && exist('aligned_str_depth_group','var') && ...
        ~isempty(aligned_str_depth_group)
        
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


