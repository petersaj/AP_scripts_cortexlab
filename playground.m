%% post to alyx

onLoad; % initializes missing-http

myAlyx = alyx.getToken([], 'andy', 'xpr1mnt');

use_date = datestr(now,'yyyy-mm-dd');
use_time = datestr(now,'HH:MM');

clear d
d.user = 'andy';

animals = {'AP030','AP031'};

water = [1.2,1.2];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.water_administered = water(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWater = alyx.postData(myAlyx, 'water-administrations', d);
end

weight = [26.7,26.3];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.weight = weight(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWeight = alyx.postData(myAlyx, 'weighings', d);
end


%% Hemidiff in V (V-V_mirror) - keep for future reference
% looks good 

U_r = reshape(U_master(:,:,1:n_vs),[],n_vs);
U_mirror_r = reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
mirror_matrix = U_r'*U_mirror_r;
fluor_allcat_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat,[],n_vs)'),size(fluor_allcat));

% sanity check: plot regular and mirrored
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_mean = squeeze(nanmean(fluor_allcat(use_trials,:,:),1))';
V_mirror_mean = squeeze(nanmean(fluor_allcat_mirror(use_trials,:,:),1))';

px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mean);
px_mirror_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mirror_mean);

AP_image_scroll([px_mean,px_mirror_mean],t);
caxis([-prctile(abs(px_mean(:)),99),prctile(abs(px_mean(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

V_hemidiff = V_mean - V_mirror_mean;
px_hemidiff = svdFrameReconstruct(U_master(:,:,1:n_vs),V_hemidiff);
AP_image_scroll(px_hemidiff,t);
caxis([-prctile(abs(px_hemidiff(:)),99),prctile(abs(px_hemidiff(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

% sanity check
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_check = squeeze(nanmean(fluor_allcat_downsamp_filt(use_trials,:,:),1))';
px_check = svdFrameReconstruct(U_master(:,:,1:n_vs),V_check);

AP_image_scroll(px_check,t);
caxis([-prctile(abs(px_check(:)),100),prctile(abs(px_check(:)),100)]);
colormap(brewermap([],'*RdBu'));
axis image;
AP_reference_outline('ccf_aligned','k');


% compatible with code: 

mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_downsamp_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_downsamp,[],n_vs)'),size(fluor_allcat_downsamp));
fluor_allcat_downsamp = fluor_allcat_downsamp - fluor_allcat_downsamp_mirror;


%% Checking waveform stuff for eventual spike sorting?

% plot templates by channel

figure; hold on;
p = arrayfun(@(x) plot(0,0,'k'),1:size(templates,3));
for curr_template = 1:size(templates,1)
    
    y = permute(templates(curr_template,:,:),[3,2,1]);
    y = y + (max(channel_positions(:,2)) - channel_positions(:,2))/1000;   
    x = (1:size(templates,2)) + channel_positions(:,1)*7;
    
    nonzero_channels = squeeze(any(templates(curr_template,:,:),2));
    [~,max_channel] = max(max(abs(templates(curr_template,:,:)),[],2),[],3);
    
    arrayfun(@(ch) set(p(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(templates,3));
    arrayfun(@(ch) set(p(ch),'Color','r'),find(nonzero_channels));
    arrayfun(@(ch) set(p(ch),'Color','k'),find(~nonzero_channels));
    set(p(max_channel),'Color','b');
    
%     ylim([min(reshape(y(nonzero_channels,:),[],1)), ...
%         max(reshape(y(nonzero_channels,:),[],1))]);
    title(curr_template);
    waitforbuttonpress;
    
end

% Get rid of zero-weight template channels
templates_permute = permute(templates,[3,2,1]);
used_template_channels = any(templates_permute,2);
if length(unique(sum(used_template_channels,1))) ~= 1
    error('Different number of unused channels');
end
templates_used = reshape( ...
    templates_permute(repmat(used_template_channels,1,size(templates_permute,2),1)), ...
    [],size(templates_permute,2),size(templates_permute,3));

template_corr = arrayfun(@(x) nanmean(AP_itril(corrcoef(templates_used(:,:,x)'),-1)),1:size(templates_used,3));
template_channel_skewness = squeeze(skewness(sum(templates_used.^2,2),[],1));
figure;plot(template_channel_skewness,template_corr,'.k');
xlabel('Channel skewness');
ylabel('Channel correlation');

figure;plot3(template_channel_skewness,template_corr,1:size(templates,1),'.k');


% 
% % Plot template skewness by depth
% skewness_time = max(skewness(templates.^2,[],2),[],3);
% skewness_channel = skewness(sum(templates.^2,2),[],3);
% 
% figure;plot3(skewness_time,skewness_channel,templateDepths,'.k')
% axis vis3d
% set(gca,'ZDir','reverse');
% xlabel('Skewness time');
% ylabel('Skewness channel');
% zlabel('Template depth');





%% Run kilosort2 on all old data

% (should be done except errored days)
% animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
% 
% for curr_animal = length(animals)
%     
%     animal = animals{curr_animal};
%     experiments = AP_find_experiments(animal);
%     ephys_days = {experiments([experiments.ephys]).day};
%     
%     for curr_day = 4
%         disp(['Kilosort2-ing: ' animal ' ' ephys_days{curr_day}])
%         AP_preprocess_phase3(animal,ephys_days{curr_day});
%     end
%       
% end

% Errored days: AP029 days 4,7,8

animal = 'AP029';
experiments = AP_find_experiments(animal);
ephys_days = {experiments([experiments.ephys]).day};

curr_day = 8;
t_range = [0,inf]; % (buffer is 65856 = 2.195, amount to leave off the end must be larger than this)
disp(['Kilosort2-ing: ' animal ' ' ephys_days{curr_day}])
AP_preprocess_phase3(animal,ephys_days{curr_day},t_range);


